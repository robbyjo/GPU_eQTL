/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.nio.file.Path;
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
    private final int maxIterations;
    private int[] cachedMask;
    private Geometry cachedGeometry;
    private int geometryBuilds;

    QPooledPedigreeModel(double[][] design,String[] groups,String pedigreeGroup,String[] individuals,
        org.jlinalg.pedigree.Pedigree pedigree,int maxIterations) {
        this.design=Arrays.stream(design).map(double[]::clone).toArray(double[][]::new);
        this.groups=groups.clone();this.individuals=individuals.clone();this.pedigreeGroup=pedigreeGroup;
        this.pedigree=pedigree;this.maxIterations=maxIterations;
        if(groups.length!=design.length||individuals.length!=design.length)throw new IllegalArgumentException("Pooled model dimensions differ");
        boolean found=false;Set<String> seen=new HashSet<>();
        for(int i=0;i<groups.length;i++) {
            if(groups[i]==null||groups[i].isBlank())throw new IllegalArgumentException("Missing residual group");
            if(groups[i].equals(pedigreeGroup)) {
                found=true;if(!seen.add(individuals[i]))throw new IllegalArgumentException("Pooled score requires one observation per pedigree individual");
                double v=pedigree.relationship(individuals[i],individuals[i]);
                if(!Double.isFinite(v)||v<=0)throw new IllegalArgumentException("Invalid pedigree diagonal");
            }
        }
        if(!found)throw new IllegalArgumentException("Pedigree group is absent from aligned participants");
        try(var algebra=SpectralLinearAlgebra.open()){engine=algebra.identity();}
    }
    static QPooledPedigreeModel create(QeQTLAnalysisConfig config,double[][] x,QCovariateTable cov,
        boolean[] selected,String[] samples)throws Exception {
        String[] ids=config.get("pedigree_id_column")==null?samples.clone():cov.columnValues(config.get("pedigree_id_column"),selected);
        if(config.get("pedigree_id_column")==null&&config.getGenotypeIdStripPrefix()!=null) {
            String prefix=config.getGenotypeIdStripPrefix();for(int i=0;i<ids.length;i++)if(ids[i].startsWith(prefix))ids[i]=ids[i].substring(prefix.length());
        }
        if(x==null){x=new double[samples.length][1];for(double[] row:x)row[0]=1;}
        return new QPooledPedigreeModel(x,cov.columnValues(config.get("residual_group_column"),selected),
            config.get("pedigree_group"),ids,QMixedModel.readPedigree(Path.of(config.getPedigreeFilename())),config.getRemlMaximumIterations());
    }
    String engine(){return engine;}
    int geometryBuilds(){return geometryBuilds;}
    String cohortCounts(int[] mask){Map<String,Integer> counts=new TreeMap<>();for(int i:mask)counts.merge(groups[i],1,Integer::sum);StringJoiner out=new StringJoiner(";");counts.forEach((name,n)->out.add(name+"="+n));return out.toString();}
    Geometry geometry(int[] mask) {
        if(!Arrays.equals(mask,cachedMask)) {
            Geometry next=new Geometry(mask);
            cachedGeometry=next;cachedMask=mask.clone();geometryBuilds++;
        }
        return cachedGeometry;
    }
    private record Rotation(int[] rows,double[] vectors) { }
    final class Geometry {
        final int[] mask;
        final int n,p;
        final double[][] q;
        final double[] eigen;
        final int[] residual;
        final List<String> names;
        final List<Rotation> rotations=new ArrayList<>();
        final double logR;
        final boolean combinedIdentity;
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
            int[] parent=new int[n];for(int i=0;i<n;i++)parent[i]=i;
            boolean nonidentity=false;
            for(int i=0;i<n;i++)if(groups[mask[i]].equals(pedigreeGroup)) {
                if(pedigree.relationship(individuals[mask[i]],individuals[mask[i]])!=1)nonidentity=true;
                for(int j=0;j<i;j++)if(groups[mask[j]].equals(pedigreeGroup)) {
                    double a=pedigree.relationship(individuals[mask[i]],individuals[mask[j]]);
                    if(!Double.isFinite(a))throw new IllegalArgumentException("Nonfinite pedigree covariance");
                    if(a!=0){parent[find(parent,i)]=find(parent,j);nonidentity=true;}
                }
            }
            combinedIdentity=!nonidentity;
            Map<Integer,List<Integer>> blocks=new LinkedHashMap<>();
            for(int i=0;i<n;i++)if(groups[mask[i]].equals(pedigreeGroup))blocks.computeIfAbsent(find(parent,i),unused->new ArrayList<>()).add(i);
            try(var algebra=SpectralLinearAlgebra.open()) {
                if(!engine.equals(algebra.identity()))throw new IllegalStateException("Pooled geometry engine changed");
                for(var block:blocks.values()) {
                    int m=block.size();int[] rows=block.stream().mapToInt(Integer::intValue).toArray();
                    if(m==1){eigen[rows[0]]=pedigree.relationship(individuals[mask[rows[0]]],individuals[mask[rows[0]]]);continue;}
                    double[] a=new double[Math.multiplyExact(m,m)];
                    for(int i=0;i<m;i++)for(int j=0;j<m;j++)a[i*m+j]=pedigree.relationship(individuals[mask[rows[i]]],individuals[mask[rows[j]]]);
                    var decomp=algebra.eigen(a,m);double scale=1;
                    for(double v:decomp.values())scale=Math.max(scale,Math.abs(v));
                    for(int i=0;i<m;i++) {
                        double value=decomp.values()[i];if(!Double.isFinite(value)||value < -1e-10*scale)throw new IllegalArgumentException("Pedigree block is not PSD");
                        eigen[rows[i]]=Math.max(0,value);
                    }
                    rotations.add(new Rotation(rows,decomp.vectors()));
                }
            }
            for(int j=0;j<p;j++){double[] column=new double[n];for(int i=0;i<n;i++)column[i]=q[i][j];column=rotate(column);for(int i=0;i<n;i++)q[i][j]=column[i];}
        }
        double[] rotate(double[] input) {
            if(input.length!=n)throw new IllegalArgumentException("Rotation dimensions differ");
            double[] result=input.clone();
            for(Rotation block:rotations){int[] rows=block.rows;int m=rows.length;
                for(int j=0;j<m;j++){double sum=0;for(int i=0;i<m;i++)sum+=block.vectors[i*m+j]*input[rows[i]];result[rows[j]]=sum;}}
            return result;
        }
        Fit fit(double[] response) {
            if(response.length!=n)throw new IllegalArgumentException("Response/mask dimensions differ");
            double[] y=rotate(response);double[] beta=new double[p];
            for(int j=0;j<p;j++)for(int i=0;i<n;i++)beta[j]+=q[i][j]*y[i];
            double sum=0;for(int i=0;i<n;i++){double v=y[i];for(int j=0;j<p;j++)v-=q[i][j]*beta[j];sum+=v*v;}
            double scale=sum/(n-p);
            if(!Double.isFinite(scale)||scale<=Math.ulp(1.0)*Math.ulp(1.0))throw new Unestimable("Constant or fixed-effect-explained response");
            double sd=Math.sqrt(scale);for(int i=0;i<n;i++)y[i]/=sd;
            int count=names.size();boolean[] active=new boolean[count];Arrays.fill(active,true);active[0]=!combinedIdentity;
            Evaluation best=null;int iterations=0;
            // Distinct deterministic initial variance allocations; no shared fit between responses.
            for(double genetic:new double[]{.1,.7}) {
                double[] theta=new double[count];Arrays.fill(theta,Math.log(.6));theta[0]=Math.log(genetic);
                if(!active[0])theta[0]=LOWER;
                Evaluation current=evaluate(theta,y,active);double[][] inverse=identity(count);
                for(int iteration=0;iteration<maxIterations;iteration++) {
                    iterations++;
                    if(current.maximumScore<=SCORE_TOLERANCE)break;
                    double[] direction=new double[count];double norm=n-p;
                    for(int a=0;a<count;a++)for(int b=0;b<count;b++)direction[a]-=inverse[a][b]*current.gradient[b]/norm;
                    for(int a=0;a<count;a++)if(!active[a]||(theta[a]<=LOWER+1e-12&&direction[a]<0)||(theta[a]>=UPPER-1e-12&&direction[a]>0))direction[a]=0;
                    double descent=dot(current.gradient,direction);
                    if(!(descent<0)){inverse=identity(count);for(int a=0;a<count;a++)direction[a]=-projected(current.gradient[a],theta[a],active[a])/norm;descent=dot(current.gradient,direction);}
                    double largest=0;for(double value:direction)largest=Math.max(largest,Math.abs(value));
                    if(largest>4){for(int a=0;a<count;a++)direction[a]*=4/largest;descent=dot(current.gradient,direction);}
                    Evaluation next=null;double[] candidate=null;
                    for(int step=0;step<45;step++) {
                        double alpha=Math.scalb(1.0,-step);candidate=theta.clone();
                        for(int a=0;a<count;a++)candidate[a]=Math.max(LOWER,Math.min(UPPER,theta[a]+alpha*direction[a]));
                        Evaluation trial=evaluate(candidate,y,active);
                        double actualDescent=0;for(int a=0;a<count;a++)actualDescent+=current.gradient[a]*(candidate[a]-theta[a]);
                        if(trial.objective<=current.objective+1e-4*actualDescent+2e-14*Math.max(1,Math.abs(current.objective))){next=trial;break;}
                    }
                    if(next==null)break;
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
                if(best==null||current.objective<best.objective)best=current;
                if(best.maximumScore<=SCORE_TOLERANCE)break;
            }
            if(best.maximumScore>SCORE_TOLERANCE)throw new ArithmeticException("Pooled REML did not converge: projected log-variance score="+best.maximumScore);
            double[] variance=best.variance.clone();for(int a=0;a<variance.length;a++)variance[a]*=scale;
            if(!active[0])variance[0]=0;
            double[] weights=best.weights.clone(),py=best.py.clone();double[][] inverse=copy(best.inverse);
            for(int i=0;i<n;i++){weights[i]/=scale;py[i]/=sd;}
            for(int a=0;a<p;a++)for(int b=0;b<p;b++)inverse[a][b]*=scale;
            double logLik=-best.objective-.5*(n-p)*Math.log(scale);
            return new Fit(this,variance,weights,py,inverse,logLik,best.maximumScore,iterations);
        }
        private Evaluation evaluate(double[] theta,double[] y,boolean[] active) {
            int count=theta.length;double[] variance=new double[count];
            for(int a=0;a<count;a++)variance[a]=active[a]?Math.exp(theta[a]):0;
            double[] w=new double[n],wy=new double[n];double[][] m=new double[p][p];double[] b=new double[p];double logdet=0;
            for(int i=0;i<n;i++) {
                double v=variance[0]*eigen[i]+variance[residual[i]];
                if(!(v>0)||!Double.isFinite(v))throw new ArithmeticException("Invalid pooled covariance");
                logdet+=Math.log(v);w[i]=1/v;wy[i]=w[i]*y[i];
                for(int a=0;a<p;a++){b[a]+=q[i][a]*wy[i];for(int c=0;c<=a;c++)m[a][c]+=q[i][a]*w[i]*q[i][c];}
            }
            for(int a=0;a<p;a++)for(int c=0;c<a;c++)m[c][a]=m[a][c];
            SmallCholesky ch=new SmallCholesky(m);double[][] inverse=ch.inverse();double[] beta=ch.solve(b),py=new double[n],gradient=new double[count];double quadratic=0;
            for(int i=0;i<n;i++) {
                double r=y[i];for(int a=0;a<p;a++)r-=q[i][a]*beta[a];py[i]=w[i]*r;quadratic+=r*py[i];
                double leverage=0;for(int a=0;a<p;a++)for(int c=0;c<p;c++)leverage+=q[i][a]*inverse[a][c]*q[i][c];
                double term=w[i]-w[i]*w[i]*leverage-py[i]*py[i];
                gradient[0]+=.5*variance[0]*eigen[i]*term;gradient[residual[i]]+=.5*variance[residual[i]]*term;
            }
            double maximum=0;for(int a=0;a<count;a++)maximum=Math.max(maximum,Math.abs(projected(gradient[a],theta[a],active[a])));
            return new Evaluation(.5*((n-p)*Math.log(2*Math.PI)+logdet+ch.logDeterminant()+2*logR+quadratic),gradient,variance,w,py,inverse,maximum);
        }
    }
    static final class Fit {
        final Geometry geometry;
        final double[] variance,weights,py;
        final double[][] inverse;
        final double logLik,maximumScore;
        final int iterations;
        Fit(Geometry geometry,double[] variance,double[] weights,double[] py,double[][] inverse,double logLik,double maximumScore,int iterations){this.geometry=geometry;this.variance=variance;this.weights=weights;this.py=py;this.inverse=inverse;this.logLik=logLik;this.maximumScore=maximumScore;this.iterations=iterations;}
        double[][] scoreColumns(){int n=geometry.n,p=geometry.p;double[][] result=new double[n][p+1];for(int i=0;i<n;i++){result[i][0]=py[i];for(int a=0;a<p;a++)result[i][a+1]=weights[i]*geometry.q[i][a];}return result;}
        double[] products(double[] rotated){double[] result=new double[geometry.p+1];for(int i=0;i<rotated.length;i++){result[0]+=rotated[i]*py[i];for(int a=0;a<geometry.p;a++)result[a+1]+=rotated[i]*weights[i]*geometry.q[i][a];}return result;}
        Association score(double[] rotated,double[] products){
            double raw=0;for(int i=0;i<rotated.length;i++)raw+=rotated[i]*rotated[i]*weights[i];
            double correction=0;for(int a=0;a<geometry.p;a++)for(int b=0;b<geometry.p;b++)correction+=products[a+1]*inverse[a][b]*products[b+1];
            double info=raw-correction;
            if(!Double.isFinite(info)||info < -1e-9*Math.max(1,raw))throw new ArithmeticException("Invalid pooled score information");
            if(info<=1e-12*Math.max(raw,1))return new Association(Double.NaN,Double.NaN,0,0,"monomorphic-or-collinear");
            double statistic=products[0]*products[0]/info;
            return new Association(products[0]/info,1/Math.sqrt(info),statistic,jdistlib.ChiSquare.cumulative(statistic,1,false,true)/Math.log(10),"ok");
        }
        String variances(){StringJoiner result=new StringJoiner(";");for(int j=0;j<variance.length;j++)result.add(geometry.names.get(j)+"="+variance[j]);return result.toString();}
    }
    record Association(double beta,double se,double statistic,double log10P,String status) { }
    private record Evaluation(double objective,double[] gradient,double[] variance,double[] weights,double[] py,double[][] inverse,double maximumScore) { }
    static final class Unestimable extends IllegalArgumentException { Unestimable(String message){super(message);} }
    static final class SmallCholesky {
        final double[][] l;
        SmallCholesky(double[][] matrix){int n=matrix.length;l=new double[n][n];for(int i=0;i<n;i++)for(int j=0;j<=i;j++){double v=matrix[i][j];for(int k=0;k<j;k++)v-=l[i][k]*l[j][k];if(i==j){if(!(v>0)||!Double.isFinite(v))throw new ArithmeticException("Pooled fixed-effect information is not positive definite");l[i][j]=Math.sqrt(v);}else l[i][j]=v/l[j][j];}}
        double[] solve(double[] b){int n=l.length;double[] x=b.clone();for(int i=0;i<n;i++){for(int j=0;j<i;j++)x[i]-=l[i][j]*x[j];x[i]/=l[i][i];}for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)x[i]-=l[j][i]*x[j];x[i]/=l[i][i];}return x;}
        double[][] inverse(){int n=l.length;double[][] inverse=new double[n][n];for(int j=0;j<n;j++){double[] b=new double[n];b[j]=1;double[] x=solve(b);for(int i=0;i<n;i++)inverse[i][j]=x[i];}return inverse;}
        double logDeterminant(){double value=0;for(int i=0;i<l.length;i++)value+=2*Math.log(l[i][i]);return value;}
    }
    private static int find(int[] parent,int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
    private static double projected(double g,double theta,boolean active){if(!active||(theta<=LOWER+1e-9&&g>0)||(theta>=UPPER-1e-9&&g<0))return 0;return g;}
    private static double dot(double[] a,double[] b){double v=0;for(int i=0;i<a.length;i++)v+=a[i]*b[i];return v;}
    private static double[][] identity(int n){double[][] v=new double[n][n];for(int i=0;i<n;i++)v[i][i]=1;return v;}
    private static double[][] copy(double[][] v){return Arrays.stream(v).map(double[]::clone).toArray(double[][]::new);}
}
