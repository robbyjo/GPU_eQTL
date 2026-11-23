/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.nio.file.Path;
import java.io.*;
import java.security.*;
import java.util.*;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.gpu.SpectralLinearAlgebra;
import gov.nih.jama.QRDecomposition;

/** Exact pedigree-block Gaussian REML. No pooled N*N covariance or projector. */
final class QPooledPedigreeModel {
    static final double SCORE_TOLERANCE=1e-6;
    private static final double LOWER=Math.log(1e-10),UPPER=Math.log(1e10);
    private final double[][] design;
    private final String[] groups,individuals;
    private final String pedigreeGroup,engine;
    private final org.jlinalg.pedigree.Pedigree pedigree;
    private final List<PedigreeBlock> pedigreeBlocks;
    private long relationshipEvaluations;
    private final int maxIterations;
    private int[] cachedMask;
    private Geometry cachedGeometry;
    private int geometryBuilds;
    private int fitExecutions;
    private String[] cacheSampleIds;

    QPooledPedigreeModel(double[][] design,String[] groups,String pedigreeGroup,String[] individuals,
        org.jlinalg.pedigree.Pedigree pedigree,int maxIterations) {
        this.design=Arrays.stream(design).map(double[]::clone).toArray(double[][]::new);
        this.groups=groups.clone();this.individuals=individuals.clone();this.pedigreeGroup=pedigreeGroup;
        cacheSampleIds=individuals.clone();
        this.pedigree=pedigree;this.maxIterations=maxIterations;
        if(groups.length!=design.length||individuals.length!=design.length)throw new IllegalArgumentException("Pooled model dimensions differ");
        boolean found=false;Set<String> seen=new HashSet<>();
        for(int i=0;i<groups.length;i++) {
            if(groups[i]==null||groups[i].isBlank())throw new IllegalArgumentException("Missing residual group");
            if(groups[i].chars().anyMatch(c->c=='\t'||c=='\r'||c=='\n'||c==';'||c=='='))throw new IllegalArgumentException("Residual group contains an audit delimiter");
            if(groups[i].equals(pedigreeGroup)) {
                found=true;if(!seen.add(individuals[i]))throw new IllegalArgumentException("Pooled score requires one observation per pedigree individual");
                double v=pedigree.relationship(individuals[i],individuals[i]);
                if(!Double.isFinite(v)||v<=0)throw new IllegalArgumentException("Invalid pedigree diagonal");
            }
        }
        if(!found)throw new IllegalArgumentException("Pedigree group is absent from aligned participants");
        try(var algebra=SpectralLinearAlgebra.open()){engine=algebra.identity();}
        pedigreeBlocks=buildPedigreeBlocks();
    }
    static QPooledPedigreeModel create(QeQTLAnalysisConfig config,double[][] x,QCovariateTable cov,
        boolean[] selected,String[] samples)throws Exception {
        String[] ids=config.get("pedigree_id_column")==null?samples.clone():cov.columnValues(config.get("pedigree_id_column"),selected);
        if(config.get("pedigree_id_column")==null&&config.getGenotypeIdStripPrefix()!=null) {
            String prefix=config.getGenotypeIdStripPrefix();for(int i=0;i<ids.length;i++)if(ids[i].startsWith(prefix))ids[i]=ids[i].substring(prefix.length());
        }
        if(x==null){x=new double[samples.length][1];for(double[] row:x)row[0]=1;}
        QPooledPedigreeModel model=new QPooledPedigreeModel(x,cov.columnValues(config.get("residual_group_column"),selected),
            config.get("pedigree_group"),ids,QMixedModel.readPedigree(Path.of(config.getPedigreeFilename())),config.getRemlMaximumIterations());
        model.cacheSampleIds=samples.clone();return model;
    }
    String engine(){return engine;}
    int geometryBuilds(){return geometryBuilds;}
    int fitExecutions(){return fitExecutions;}
    long relationshipEvaluations(){return relationshipEvaluations;}
    int pedigreeBlockCount(){return pedigreeBlocks.size();}
    int largestPedigreeBlock(){return pedigreeBlocks.stream().mapToInt(block->block.samples.length).max().orElse(0);}
    String cohortCounts(int[] mask){Map<String,Integer> counts=new TreeMap<>();for(int i:mask)counts.merge(groups[i],1,Integer::sum);StringJoiner out=new StringJoiner(";");counts.forEach((name,n)->out.add(name+"="+n));return out.toString();}
    Geometry geometry(int[] mask) {
        if(!Arrays.equals(mask,cachedMask)) {
            Geometry next=new Geometry(mask);
            cachedGeometry=next;cachedMask=mask.clone();geometryBuilds++;
        }
        return cachedGeometry;
    }
    private record Rotation(int[] rows,double[] vectors) { }
    private record PedigreeBlock(int[] samples,double[] covariance) { }
    private double relationship(int a,int b){relationshipEvaluations++;return pedigree.relationship(individuals[a],individuals[b]);}
    private List<PedigreeBlock> buildPedigreeBlocks() {
        int[] selected=java.util.stream.IntStream.range(0,groups.length).filter(i->groups[i].equals(pedigreeGroup)).toArray();
        int[] parent=new int[selected.length];for(int i=0;i<parent.length;i++)parent[i]=i;
        for(int i=0;i<selected.length;i++)for(int j=0;j<i;j++) {
            double value=relationship(selected[i],selected[j]);
            if(!Double.isFinite(value))throw new IllegalArgumentException("Nonfinite pedigree covariance");
            if(value!=0)parent[find(parent,i)]=find(parent,j);
        }
        Map<Integer,List<Integer>> components=new LinkedHashMap<>();
        for(int i=0;i<selected.length;i++)components.computeIfAbsent(find(parent,i),unused->new ArrayList<>()).add(selected[i]);
        List<PedigreeBlock> result=new ArrayList<>();long bytes=0;
        for(List<Integer> component:components.values()) {
            int[] samples=component.stream().mapToInt(Integer::intValue).toArray();int m=samples.length;
            bytes=Math.addExact(bytes,Math.multiplyExact(Math.multiplyExact((long)m,m),8));
            if(bytes>Runtime.getRuntime().maxMemory()*.3)throw new IllegalArgumentException("Pedigree component cache exceeds 30% of JVM heap; increase -Xmx");
            double[] covariance=new double[Math.multiplyExact(m,m)];
            for(int i=0;i<m;i++)for(int j=0;j<=i;j++) {
                double value=relationship(samples[i],samples[j]);
                if(!Double.isFinite(value))throw new IllegalArgumentException("Nonfinite pedigree covariance");
                covariance[i*m+j]=covariance[j*m+i]=value;
            }
            result.add(new PedigreeBlock(samples,covariance));
        }
        return List.copyOf(result);
    }
    final class Geometry {
        final int[] mask;
        final int n,p;
        final double[][] q;
        final double[] qFlat;
        final double[] eigen;
        final int[] residual;
        final List<String> names;
        final List<Rotation> rotations=new ArrayList<>();
        final double logR;
        final boolean combinedIdentity;
        final double[] lowerBounds,upperBounds;
        Geometry(int[] mask) {
            this.mask=mask.clone();n=mask.length;
            Set<Integer> unique=new HashSet<>();double[][] local=new double[n][];
            TreeSet<String> levels=new TreeSet<>();
            for(int i=0;i<n;i++) {
                if(mask[i]<0||mask[i]>=design.length||!unique.add(mask[i]))throw new IllegalArgumentException("Invalid pooled sample mask");
                local[i]=design[mask[i]].clone();levels.add(groups[mask[i]]);
            }
            if(n<3)throw new Unestimable("Too few observed samples");
            QRDecomposition qr=new QRDecomposition(local);int rank=qr.getRank();
            if(n<=rank+1)throw new Unestimable("No positive association residual DF");
            // Preserve the exact column space if a mask loses a nuisance-factor level.
            if(rank<local[0].length) {
                int[] pivot=qr.getPivot();double[][] reduced=new double[n][rank];
                for(int i=0;i<n;i++)for(int j=0;j<rank;j++)reduced[i][j]=local[i][pivot[j]];
                qr=new QRDecomposition(reduced);
            }
            p=rank;q=qr.getQArray();double[][] r=qr.getRArray();double log=0;
            for(int j=0;j<p;j++)log+=Math.log(Math.abs(r[j][j]));logR=log;
            eigen=new double[n];residual=new int[n];
            List<String> sorted=new ArrayList<>(levels);names=new ArrayList<>();names.add("additive:"+pedigreeGroup);
            for(String level:sorted)names.add("residual:"+level);
            for(int i=0;i<n;i++)residual[i]=1+sorted.indexOf(groups[mask[i]]);
            int[] localIndex=new int[groups.length];Arrays.fill(localIndex,-1);
            for(int i=0;i<n;i++)localIndex[mask[i]]=i;
            boolean proportional=true;double diagonal=Double.NaN;
            try(var algebra=SpectralLinearAlgebra.open()) {
                if(!engine.equals(algebra.identity()))throw new IllegalStateException("Pooled geometry engine changed");
                for(PedigreeBlock block:pedigreeBlocks) {
                    int full=block.samples.length;
                    int[] selected=java.util.stream.IntStream.range(0,full).filter(i->localIndex[block.samples[i]]>=0).toArray();
                    int m=selected.length;if(m==0)continue;
                    int[] rows=Arrays.stream(selected).map(i->localIndex[block.samples[i]]).toArray();
                    for(int i=0;i<m;i++) {
                        double value=block.covariance[selected[i]*full+selected[i]];
                        if(Double.isNaN(diagonal))diagonal=value;else if(value!=diagonal)proportional=false;
                        for(int j=0;j<i;j++)if(block.covariance[selected[i]*full+selected[j]]!=0)proportional=false;
                    }
                    if(m==1){eigen[rows[0]]=block.covariance[selected[0]*full+selected[0]];continue;}
                    double[] a=new double[Math.multiplyExact(m,m)];
                    for(int i=0;i<m;i++)for(int j=0;j<m;j++)a[i*m+j]=block.covariance[selected[i]*full+selected[j]];
                    var decomp=algebra.eigen(a,m);double scale=1;
                    for(double v:decomp.values())scale=Math.max(scale,Math.abs(v));
                    for(int i=0;i<m;i++) {
                        double value=decomp.values()[i];if(!Double.isFinite(value)||value < -1e-10*scale)throw new IllegalArgumentException("Pedigree block is not PSD");
                        eigen[rows[i]]=Math.max(0,value);
                    }
                    rotations.add(new Rotation(rows,decomp.vectors()));
                }
            }
            combinedIdentity=proportional;
            lowerBounds=new double[names.size()];upperBounds=new double[names.size()];Arrays.fill(lowerBounds,LOWER);Arrays.fill(upperBounds,UPPER);
            if(combinedIdentity) {
                names.set(0,"unidentified-additive:"+pedigreeGroup);
                if(!Double.isNaN(diagonal)) {
                    int index=1+sorted.indexOf(pedigreeGroup);double shift=Math.log1p(diagonal);
                    lowerBounds[index]+=shift;upperBounds[index]+=shift;names.set(index,"total:"+pedigreeGroup);
                }
            }
            for(int j=0;j<p;j++){double[] column=new double[n];for(int i=0;i<n;i++)column[i]=q[i][j];column=rotate(column);for(int i=0;i<n;i++)q[i][j]=column[i];}
            qFlat=new double[Math.multiplyExact(n,p)];for(int i=0;i<n;i++)System.arraycopy(q[i],0,qFlat,i*p,p);
            requireIdentifiableCovariance();
        }
        private void requireIdentifiableCovariance() {
            int count=names.size();double[][][] projected=new double[count][p][p];double[][] gram=new double[count][count];double[] rawNorm=new double[count];
            for(int i=0;i<n;i++) {
                double leverage=dot(q[i],q[i]);double[] kernel=new double[count];kernel[0]=eigen[i];kernel[residual[i]]=1;
                for(int a=0;a<count;a++)if(kernel[a]!=0) {
                    rawNorm[a]+=kernel[a]*kernel[a];
                    for(int b=0;b<count;b++)gram[a][b]+=kernel[a]*kernel[b]*(1-2*leverage);
                    for(int j=0;j<p;j++)for(int k=0;k<p;k++)projected[a][j][k]+=kernel[a]*q[i][j]*q[i][k];
                }
            }
            for(int a=0;a<count;a++)for(int b=0;b<count;b++)for(int j=0;j<p;j++)for(int k=0;k<p;k++)gram[a][b]+=projected[a][j][k]*projected[b][j][k];
            int first=combinedIdentity?1:0;double[][] normalized=new double[count-first][count-first];
            for(int a=first;a<count;a++) {
                if(gram[a][a]<=1e-10*Math.max(1,rawNorm[a]))throw new Unestimable("Covariance component has no estimable residual contrasts: "+names.get(a));
                for(int b=first;b<count;b++)normalized[a-first][b-first]=gram[a][b]/Math.sqrt(gram[a][a]*gram[b][b]);
            }
            try {
                SmallCholesky factor=new SmallCholesky(normalized);
                for(int i=0;i<factor.l.length;i++)if(factor.l[i][i]*factor.l[i][i]<=1e-10)throw new Unestimable("Covariance components are not identifiable after fixed-effect projection");
            }catch(ArithmeticException error){throw new Unestimable("Covariance components are not identifiable after fixed-effect projection");}
        }
        double[] rotate(double[] input) {
            if(input.length!=n)throw new IllegalArgumentException("Rotation dimensions differ");
            double[] result=input.clone();
            for(Rotation block:rotations){int[] rows=block.rows;int m=rows.length;
                for(int j=0;j<m;j++){double sum=0;for(int i=0;i<m;i++)sum+=block.vectors[i*m+j]*input[rows[i]];result[rows[j]]=sum;}}
            return result;
        }
        double[][] rotateRows(double[][] inputs) {
            double[][] result=Arrays.stream(inputs).map(double[]::clone).toArray(double[][]::new);
            for(double[] input:inputs)if(input.length!=n)throw new IllegalArgumentException("Genotype rotation dimensions differ");
            try(var algebra=SpectralLinearAlgebra.open()) {
                if(!engine.equals(algebra.identity()))throw new IllegalStateException("Pooled rotation engine changed");
                for(Rotation block:rotations) {
                    int m=block.rows.length;
                    if(m<16) {
                        for(int a=0;a<inputs.length;a++)for(int j=0;j<m;j++){double value=0;for(int i=0;i<m;i++)value+=inputs[a][block.rows[i]]*block.vectors[i*m+j];result[a][block.rows[j]]=value;}
                    } else {
                        double[] selected=new double[Math.multiplyExact(inputs.length,m)];
                        for(int a=0;a<inputs.length;a++)for(int i=0;i<m;i++)selected[a*m+i]=inputs[a][block.rows[i]];
                        double[] rotated=algebra.multiply(selected,false,block.vectors,false,inputs.length,m,m);
                        for(int a=0;a<inputs.length;a++)for(int j=0;j<m;j++)result[a][block.rows[j]]=rotated[a*m+j];
                    }
                }
            }
            return result;
        }
        String nullCacheKey(String trait,double[] response)throws IOException {
            try {
                MessageDigest digest=MessageDigest.getInstance("SHA-256");
                try(DataOutputStream out=new DataOutputStream(new DigestOutputStream(OutputStream.nullOutputStream(),digest))) {
                    out.writeUTF("pooled-null-v5-relative-reml");out.writeUTF(engine);out.writeUTF(trait);
                    out.writeInt(maxIterations);out.writeInt(n);out.writeInt(p);out.writeDouble(logR);out.writeBoolean(combinedIdentity);
                    for(int i:mask){out.writeInt(i);out.writeUTF(cacheSampleIds[i]);out.writeUTF(groups[i]);out.writeUTF(individuals[i]);}
                    for(String name:names)out.writeUTF(name);
                    for(double value:eigen)out.writeDouble(value);for(int value:residual)out.writeInt(value);
                    for(double value:lowerBounds)out.writeDouble(value);for(double value:upperBounds)out.writeDouble(value);
                    for(double[] row:q)for(double value:row)out.writeDouble(value);
                    for(Rotation rotation:rotations){out.writeInt(rotation.rows.length);for(int row:rotation.rows)out.writeInt(row);for(double value:rotation.vectors)out.writeDouble(value);}
                    for(double value:response)out.writeDouble(value);
                }
                return HexFormat.of().formatHex(digest.digest());
            }catch(NoSuchAlgorithmException impossible){throw new IllegalStateException(impossible);}
        }
        Fit fit(double[] response) {
            fitExecutions++;
            try(var algebra=SpectralLinearAlgebra.open()) {
                if(!engine.equals(algebra.identity()))throw new IllegalStateException("Pooled fitter engine changed");
                return fitWithAlgebra(response,algebra);
            }
        }
        private Fit fitWithAlgebra(double[] response,SpectralLinearAlgebra algebra) {
            if(response.length!=n)throw new IllegalArgumentException("Response/mask dimensions differ");
            double magnitude=0;for(double value:response){if(!Double.isFinite(value))throw new IllegalArgumentException("Nonfinite response");magnitude=Math.max(magnitude,Math.abs(value));}
            if(magnitude==0)throw new Unestimable("Constant or fixed-effect-explained response");
            double[] normalized=response.clone();for(int i=0;i<n;i++)normalized[i]/=magnitude;
            double[] y=rotate(normalized);double[] beta=new double[p];
            for(int j=0;j<p;j++)for(int i=0;i<n;i++)beta[j]+=q[i][j]*y[i];
            double sum=0,maximumResponse=0;for(int i=0;i<n;i++){if(!Double.isFinite(y[i]))throw new IllegalArgumentException("Nonfinite pooled response");maximumResponse=Math.max(maximumResponse,Math.abs(y[i]));double v=y[i];for(int j=0;j<p;j++)v-=q[i][j]*beta[j];y[i]=v;sum+=v*v;}
            double scale=sum/(n-p);
            if(!Double.isFinite(scale)||scale<=Math.ulp(1.0)*Math.ulp(1.0))throw new Unestimable("Constant or fixed-effect-explained response");
            double sd=Math.sqrt(scale);for(int i=0;i<n;i++)y[i]/=sd;
            double responseScale=magnitude*sd;
            int count=names.size();boolean[] active=new boolean[count];Arrays.fill(active,true);active[0]=!combinedIdentity;
            Evaluation best=null;int iterations=0;
            // Distinct deterministic initial variance allocations; no shared fit between responses.
            int startNumber=0;
            for(double genetic:new double[]{.1,.7}) {
                int budget=startNumber++==0?Math.max(1,maxIterations/2):maxIterations-iterations;
                double[] theta=new double[count];Arrays.fill(theta,Math.log(.6));theta[0]=Math.log(genetic);
                if(!active[0])theta[0]=LOWER;
                Evaluation current=evaluate(theta,y,active,algebra);double[][] inverse=identity(count);
                for(int iteration=0;iteration<budget;iteration++) {
                    iterations++;
                    if(current.maximumScore<=SCORE_TOLERANCE)break;
                    double[] direction=new double[count];double norm=n-p;
                    for(int a=0;a<count;a++)for(int b=0;b<count;b++)direction[a]-=inverse[a][b]*current.gradient[b]/norm;
                    for(int a=0;a<count;a++)if(!active[a]||(theta[a]<=lowerBounds[a]+1e-12&&direction[a]<0)||(theta[a]>=upperBounds[a]-1e-12&&direction[a]>0))direction[a]=0;
                    double descent=dot(current.gradient,direction);
                    if(!(descent<0)){inverse=identity(count);for(int a=0;a<count;a++)direction[a]=-projected(current.gradient[a],theta[a],active[a],lowerBounds[a],upperBounds[a])/norm;descent=dot(current.gradient,direction);}
                    double largest=0;for(double value:direction)largest=Math.max(largest,Math.abs(value));
                    if(largest>4){for(int a=0;a<count;a++)direction[a]*=4/largest;descent=dot(current.gradient,direction);}
                    Evaluation next=null;double[] candidate=null;
                    for(int step=0;step<45;step++) {
                        double alpha=Math.scalb(1.0,-step);candidate=theta.clone();
                        for(int a=0;a<count;a++)candidate[a]=Math.max(lowerBounds[a],Math.min(upperBounds[a],theta[a]+alpha*direction[a]));
                        Evaluation trial=evaluate(candidate,y,active,algebra);
                        double actualDescent=0;for(int a=0;a<count;a++)actualDescent+=current.gradient[a]*(candidate[a]-theta[a]);
                        double change=objectiveChange(current,theta,candidate,algebra);
                        if(change<=1e-4*actualDescent+2e-14*Math.max(1,Math.abs(current.objective))){next=trial;break;}
                    }
                    if(next==null||Arrays.equals(theta,candidate))break;
                    double[] s=new double[count],z=new double[count];
                    for(int a=0;a<count;a++){s[a]=candidate[a]-theta[a];z[a]=(next.gradient[a]-current.gradient[a])/norm;}
                    double sz=dot(s,z);
                    if(sz>1e-14*dot(s,s)) {
                        double[] hz=new double[count];for(int a=0;a<count;a++)for(int b=0;b<count;b++)hz[a]+=inverse[a][b]*z[b];
                        double coefficient=(sz+dot(z,hz))/(sz*sz);
                        for(int a=0;a<count;a++)for(int b=0;b<count;b++)inverse[a][b]+=coefficient*s[a]*s[b]-(hz[a]*s[b]+s[a]*hz[b])/sz;
                    } else inverse=identity(count);
                    theta=candidate;current=next;
                }
                // A stationary restart is valid even if roundoff makes an earlier,
                // nonstationary objective marginally smaller. Never discard it.
                if(current.maximumScore<=SCORE_TOLERANCE){best=current;break;}
                if(best==null||current.objective<best.objective)best=current;
            }
            if(!Double.isFinite(best.maximumScore)||best.maximumScore>SCORE_TOLERANCE)throw new ArithmeticException("Pooled REML did not converge: projected log-variance score="+best.maximumScore);
            double[] variance=best.variance.clone();for(int a=0;a<variance.length;a++){variance[a]=(variance[a]*responseScale)*responseScale;if(!Double.isFinite(variance[a])||(active[a]&&variance[a]==0))throw new ArithmeticException("Pooled variance exceeds FP64 range");}
            if(!active[0])variance[0]=0;
            double[] weights=best.weights.clone(),py=best.py.clone();double[][] inverse=copy(best.inverse);
            // Preserve v1.9.2 normalized score state; convert only reported effects/SE.
            double logLik=-best.objective-(n-p)*Math.log(responseScale);
            return new Fit(this,variance,weights,py,inverse,responseScale,logLik,best.maximumScore,iterations);
        }
        Objective objective(double[] logVariances,double[] response) {
            if(logVariances.length!=names.size()||response.length!=n)throw new IllegalArgumentException("Diagnostic dimensions differ");
            boolean[] active=new boolean[names.size()];Arrays.fill(active,true);active[0]=!combinedIdentity;
            try(var algebra=SpectralLinearAlgebra.open()) {
                Evaluation result=evaluate(logVariances,rotate(response),active,algebra);
                return new Objective(result.objective,result.gradient.clone());
            }
        }
        double objectiveChange(double[] from,double[] to,double[] response) {
            if(from.length!=names.size()||to.length!=names.size()||response.length!=n)throw new IllegalArgumentException("Diagnostic dimensions differ");
            boolean[] active=new boolean[names.size()];Arrays.fill(active,true);active[0]=!combinedIdentity;
            try(var algebra=SpectralLinearAlgebra.open()) {
                return objectiveChange(evaluate(from,rotate(response),active,algebra),from,to,algebra);
            }
        }
        // Evaluate a likelihood difference directly, without subtracting two
        // independently rounded weighted-QR objectives near a stationary point.
        private double objectiveChange(Evaluation current,double[] theta,double[] candidate,SpectralLinearAlgebra algebra) {
            int count=theta.length;double[] delta=new double[count];
            for(int a=0;a<count;a++)delta[a]=current.variance[a]*Math.expm1(candidate[a]-theta[a]);
            double[] shifted=new double[Math.multiplyExact(n,p)],cross=new double[p];double covariance=0,quadratic=0;
            for(int i=0;i<n;i++) {
                double ratio=current.weights[i]*(delta[0]*eigen[i]+delta[residual[i]]);
                if(!(ratio>-1)||!Double.isFinite(ratio))throw new ArithmeticException("Invalid relative pooled covariance");
                double change=-ratio/(1+ratio),r=current.py[i]/Math.sqrt(current.weights[i]);
                covariance+=Math.log1p(ratio);quadratic+=change*r*r;
                for(int a=0;a<p;a++){double value=change*current.orthogonal[i*p+a];shifted[i*p+a]=value;cross[a]+=value*r;}
            }
            double[] flat=algebra.multiply(current.orthogonal,true,shifted,false,p,p,n);
            double[][] increment=new double[p][p];for(int a=0;a<p;a++)System.arraycopy(flat,a*p,increment[a],0,p);
            // With A=orth(sqrt(W)X), T=W_new/W-1 and e=sqrt(W)r:
            // delta SSE=e'Te-(A'Te)'(I+A'TA)^-1(A'Te).
            SmallCholesky relative=new SmallCholesky(increment,true);
            double value=.5*(covariance+relative.logDeterminant()+quadratic-dot(cross,relative.solve(cross)));
            if(!Double.isFinite(value))throw new ArithmeticException("Nonfinite pooled REML objective change");
            return value;
        }
        private Evaluation evaluate(double[] theta,double[] y,boolean[] active,SpectralLinearAlgebra algebra) {
            int count=theta.length;double[] variance=new double[count];
            for(int a=0;a<count;a++)variance[a]=active[a]?Math.exp(theta[a]):0;
            double[] w=new double[n],sqrtW=new double[n],wy=new double[n],weightedQ=new double[Math.multiplyExact(n,p)];double logdet=0;
            for(int i=0;i<n;i++) {
                double v=variance[0]*eigen[i]+variance[residual[i]];
                if(!(v>0)||!Double.isFinite(v))throw new ArithmeticException("Invalid pooled covariance");
                logdet+=Math.log(v);w[i]=1/v;sqrtW[i]=Math.sqrt(w[i]);wy[i]=sqrtW[i]*y[i];
                for(int a=0;a<p;a++)weightedQ[i*p+a]=sqrtW[i]*q[i][a];
            }
            // Whiten, then project through thin Householder QR. Forming Q'WQ
            // squares the condition number when one cohort reaches a variance bound.
            var projection=algebra.projectResidual(weightedQ,wy,n,p);
            double[] upper=projection.upper(),inverseR=new double[p*p];double logInformation=0;
            for(int a=0;a<p;a++)logInformation+=2*Math.log(Math.abs(upper[a*p+a]));
            for(int col=0;col<p;col++)for(int row=p-1;row>=0;row--) {
                double value=row==col?1:0;
                for(int k=row+1;k<p;k++)value-=upper[row*p+k]*inverseR[k*p+col];
                inverseR[row*p+col]=value/upper[row*p+row];
            }
            double[] inverseFlat=algebra.multiply(inverseR,false,inverseR,true,p,p,p);
            double[][] inverse=new double[p][p];for(int a=0;a<p;a++)System.arraycopy(inverseFlat,a*p,inverse[a],0,p);
            double[] py=new double[n],gradient=new double[count];double quadratic=0;
            for(int i=0;i<n;i++) {
                double r=projection.residual()[i];py[i]=sqrtW[i]*r;quadratic+=r*r;
                double leverage=0;for(int a=0;a<p;a++){double value=projection.orthonormalColumns()[i*p+a];leverage+=value*value;}
                double term=w[i]*(1-leverage)-py[i]*py[i];
                gradient[0]+=.5*variance[0]*eigen[i]*term;gradient[residual[i]]+=.5*variance[residual[i]]*term;
            }
            double maximum=0;for(int a=0;a<count;a++)maximum=Math.max(maximum,Math.abs(projected(gradient[a],theta[a],active[a],lowerBounds[a],upperBounds[a])));
            double objective=.5*((n-p)*Math.log(2*Math.PI)+logdet+logInformation+2*logR+quadratic);
            if(!Double.isFinite(objective)||!Double.isFinite(maximum))throw new ArithmeticException("Nonfinite pooled REML objective or score");
            return new Evaluation(objective,gradient,variance,w,py,inverse,maximum,projection.orthonormalColumns());
        }
    }
    static final class Fit {
        final Geometry geometry;
        final double[] variance,weights,py;
        final double[][] inverse;
        private final double[] sqrtWeights,scoreBasis;
        final double responseScale,logLik,maximumScore;
        final int iterations;
        final boolean reused;
        Fit(Geometry geometry,double[] variance,double[] weights,double[] py,double[][] inverse,double responseScale,double logLik,double maximumScore,int iterations){this(geometry,variance,weights,py,inverse,responseScale,logLik,maximumScore,iterations,false);}
        Fit(Geometry geometry,double[] variance,double[] weights,double[] py,double[][] inverse,double responseScale,double logLik,double maximumScore,int iterations,boolean reused){
            this.geometry=geometry;this.variance=variance;this.weights=weights;this.py=py;this.inverse=inverse;this.responseScale=responseScale;this.logLik=logLik;this.maximumScore=maximumScore;this.iterations=iterations;this.reused=reused;
            int n=geometry.n,p=geometry.p;sqrtWeights=new double[n];double[] design=new double[Math.multiplyExact(n,p)];
            for(int i=0;i<n;i++){sqrtWeights[i]=Math.sqrt(weights[i]);for(int a=0;a<p;a++)design[i*p+a]=sqrtWeights[i]*geometry.q[i][a];}
            // Reconstruct from the immutable null. The
            // fitted model and serialized state do not change with this score basis.
            try(var algebra=SpectralLinearAlgebra.open()){scoreBasis=algebra.projectResidual(design,new double[n],n,p).orthonormalColumns();}
        }
        double[][] scoreColumns(){int n=geometry.n,p=geometry.p;double[][] result=new double[n][p+1];for(int i=0;i<n;i++){result[i][0]=py[i];for(int a=0;a<p;a++)result[i][a+1]=sqrtWeights[i]*scoreBasis[i*p+a];}return result;}
        double[] products(double[] rotated){int p=geometry.p;double[] result=new double[p+1];for(int i=0;i<rotated.length;i++){result[0]+=rotated[i]*py[i];for(int a=0;a<p;a++)result[a+1]+=rotated[i]*sqrtWeights[i]*scoreBasis[i*p+a];}return result;}
        Association score(double[] rotated,double[] products){
            double raw=0;for(int i=0;i<rotated.length;i++){double value=rotated[i]*sqrtWeights[i];raw+=value*value;}
            double correction=0;for(int a=0;a<geometry.p;a++)correction+=products[a+1]*products[a+1];
            double info=raw-correction,score=products[0];
            if(!Double.isFinite(info)||info < -1e-9*raw)throw new ArithmeticException("Invalid pooled score information");
            if(info<=.01*raw&&raw>0) {
                // Close subtraction can lose the small positive information at a
                // cohort variance bound. Use its equivalent residual norm and
                // reproject once to remove roundoff along the fixed-effect basis.
                int n=geometry.n,p=geometry.p;double[] residual=new double[n],coefficients=new double[p];
                for(int i=0;i<n;i++){double value=rotated[i]*sqrtWeights[i];for(int a=0;a<p;a++)value-=scoreBasis[i*p+a]*products[a+1];residual[i]=value;}
                for(int i=0;i<n;i++)for(int a=0;a<p;a++)coefficients[a]+=scoreBasis[i*p+a]*residual[i];
                info=0;score=0;
                for(int i=0;i<n;i++){double value=residual[i];for(int a=0;a<p;a++)value-=scoreBasis[i*p+a]*coefficients[a];info+=value*value;score+=value*(py[i]/sqrtWeights[i]);}
            }
            if(info<=1e-12*raw)return new Association(Double.NaN,Double.NaN,0,0,"monomorphic-or-collinear");
            double z=score/Math.sqrt(info),statistic=z*z;
            if(!Double.isFinite(statistic))throw new ArithmeticException("Nonfinite pooled score statistic");
            double beta=(score/info)*responseScale,se=responseScale/Math.sqrt(info);
            if(!Double.isFinite(beta)||!Double.isFinite(se)||se<=0)throw new ArithmeticException("Pooled effect or SE exceeds FP64 range");
            return new Association(beta,se,statistic,jdistlib.ChiSquare.cumulative(statistic,1,false,true)/Math.log(10),"ok");
        }
        String variances(){StringJoiner result=new StringJoiner(";");for(int j=0;j<variance.length;j++)result.add(geometry.names.get(j)+"="+variance[j]);return result.toString();}
    }
    record Association(double beta,double se,double statistic,double log10P,String status) { }
    record Objective(double value,double[] gradient) { }
    private record Evaluation(double objective,double[] gradient,double[] variance,double[] weights,double[] py,double[][] inverse,double maximumScore,double[] orthogonal) { }
    static final class Unestimable extends IllegalArgumentException { Unestimable(String message){super(message);} }
    static final class SmallCholesky {
        final double[][] l;private final double logDeterminant;
        SmallCholesky(double[][] matrix){this(matrix,false);}
        SmallCholesky(double[][] matrix,boolean identityPlus){
            int n=matrix.length;l=new double[n][n];double determinant=0;
            for(int i=0;i<n;i++)for(int j=0;j<=i;j++) {
                double v=matrix[i][j];for(int k=0;k<j;k++)v-=l[i][k]*l[j][k];
                if(i==j){double diagonal=identityPlus?1+v:v;if(!(diagonal>0)||!Double.isFinite(diagonal))throw new ArithmeticException("Pooled fixed-effect information is not positive definite");l[i][j]=Math.sqrt(diagonal);determinant+=identityPlus?Math.log1p(v):Math.log(v);}
                else l[i][j]=v/l[j][j];
            }
            logDeterminant=determinant;
        }
        double[] solve(double[] b){int n=l.length;double[] x=b.clone();for(int i=0;i<n;i++){for(int j=0;j<i;j++)x[i]-=l[i][j]*x[j];x[i]/=l[i][i];}for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)x[i]-=l[j][i]*x[j];x[i]/=l[i][i];}return x;}
        double[][] inverse(){int n=l.length;double[][] inverse=new double[n][n];for(int j=0;j<n;j++){double[] b=new double[n];b[j]=1;double[] x=solve(b);for(int i=0;i<n;i++)inverse[i][j]=x[i];}return inverse;}
        double logDeterminant(){return logDeterminant;}
    }
    private static int find(int[] parent,int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
    private static double projected(double g,double theta,boolean active,double lower,double upper){if(!active||(theta<=lower+1e-9&&g>0)||(theta>=upper-1e-9&&g<0))return 0;return g;}
    private static double dot(double[] a,double[] b){double v=0;for(int i=0;i<a.length;i++)v+=a[i]*b[i];return v;}
    private static double[][] identity(int n){double[][] v=new double[n][n];for(int i=0;i<n;i++)v[i][i]=1;return v;}
    private static double[][] copy(double[][] v){return Arrays.stream(v).map(double[]::clone).toArray(double[][]::new);}
}
