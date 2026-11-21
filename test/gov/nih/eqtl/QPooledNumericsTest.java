/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledNumericsTest {
    @Test void objectiveAndAnalyticScoresMatchIndependentDenseRAndFiniteDifferences()throws Exception {
        var d=QPooledPedigreeModelTest.data();var model=QPooledPedigreeModelTest.model(d);
        for(String line:Files.readAllLines(QPooledPedigreeModelTest.ROOT.resolve("diagnostics.csv")).subList(1,3)) {
            String[] expected=line.split(",");Set<Integer> missing=expected[0].equals("missing")?Set.of(2,7,18,41,52,80,93):Set.of();
            int[] mask=java.util.stream.IntStream.range(0,96).filter(i->!missing.contains(i)).toArray();
            var geometry=model.geometry(mask);double[] y=Arrays.stream(mask).mapToDouble(i->d.y()[i]).toArray();
            double[] point=Arrays.stream(new double[]{.2,.03,.15,.34}).map(Math::log).toArray();var actual=geometry.objective(point,y);
            assertEquals(Double.parseDouble(expected[1]),actual.value(),1e-9);
            for(int j=0;j<4;j++) {
                assertEquals(Double.parseDouble(expected[j+2]),actual.gradient()[j],2e-9);
                double[] plus=point.clone(),minus=point.clone();plus[j]+=1e-5;minus[j]-=1e-5;
                assertEquals(actual.gradient()[j],(geometry.objective(plus,y).value()-geometry.objective(minus,y).value())/2e-5,1e-7);
            }
        }
    }
    @Test void nativeFitMatchesDenseR()throws Exception {
        String previous=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","openblas");
        try {
            var d=QPooledPedigreeModelTest.data();var geometry=QPooledPedigreeModelTest.model(d).geometry(java.util.stream.IntStream.range(0,96).toArray());var fit=geometry.fit(d.y());
            String[] expected=Files.readAllLines(QPooledPedigreeModelTest.ROOT.resolve("expected.csv")).get(1).split(",");double[] g=geometry.rotate(d.g());var score=fit.score(g,fit.products(g));
            assertEquals(Double.parseDouble(expected[10]),fit.logLik,2e-5);assertEquals(Double.parseDouble(expected[9]),score.log10P(),2e-6);
            assertTrue(fit.maximumScore<=1e-6);
        }finally{if(previous==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",previous);}
    }
    @Test void confoundedCovariancesAndNoContrastCohortsAreExplicitlyUnestimable()throws Exception {
        var d=QPooledPedigreeModelTest.data();double[][] x=new double[96][15];
        for(int i=0;i<96;i++){x[i][0]=1;x[i][1]=i>=48&&i<72?1:0;x[i][2]=i>=72?1:0;if(i<48)x[i][3+i/4]=1;}
        var model=new QPooledPedigreeModel(x,d.group(),"FHS",d.ids(),QMixedModel.readPedigree(QPooledPedigreeModelTest.ROOT.resolve("pedigree.tsv")),1000);
        assertThrows(QPooledPedigreeModel.Unestimable.class,()->model.geometry(java.util.stream.IntStream.range(0,96).toArray()));
        var original=QPooledPedigreeModelTest.model(d);int[] oneCardia=java.util.stream.IntStream.range(0,73).toArray();
        assertThrows(QPooledPedigreeModel.Unestimable.class,()->original.geometry(oneCardia));
    }
    @Test void proportionalInbredPedigreeRetainsSummedVarianceBounds(@TempDir Path root)throws Exception {
        StringBuilder ped=new StringBuilder("ID\tSIRE\tDAM\n");
        for(int f=0;f<2;f++)ped.append("f"+f+"\t0\t0\nm"+f+"\t0\t0\na"+f+"\tf"+f+"\tm"+f+"\nb"+f+"\tf"+f+"\tm"+f+"\nc"+f+"\ta"+f+"\tb"+f+"\n");
        Path file=root.resolve("ped.tsv");Files.writeString(file,ped);var pedigree=QMixedModel.readPedigree(file);assertEquals(1.25,pedigree.relationship("c0","c0"),0);
        String[] ids={"c0","c1","u0","u1","u2","u3","u4","u5"},groups={"FHS","FHS","WHI","WHI","WHI","WHI","WHI","WHI"};
        double[][] x=new double[8][2];for(int i=0;i<8;i++){x[i][0]=1;x[i][1]=i>=2?1:0;}
        var geometry=new QPooledPedigreeModel(x,groups,"FHS",ids,pedigree,1000).geometry(java.util.stream.IntStream.range(0,8).toArray());
        assertTrue(geometry.combinedIdentity);assertEquals(Math.log(2.25e-10),geometry.lowerBounds[1],1e-14);assertEquals(Math.log(2.25e10),geometry.upperBounds[1],1e-14);
        double[] y={1,3,2,4,5,6,8,9};var fit=geometry.fit(y);double mean=Arrays.stream(y,2,8).average().orElseThrow(),variance=0;
        for(int i=2;i<8;i++)variance+=(y[i]-mean)*(y[i]-mean)/5;
        assertEquals(2,fit.variance[1],2e-5);assertEquals(variance,fit.variance[2],2e-5);assertEquals(0,fit.variance[0]);
    }
    @Test void nullCacheReusesAcrossModelsAndRejectsCorruption(@TempDir Path root)throws Exception {
        var d=QPooledPedigreeModelTest.data();int[] mask=java.util.stream.IntStream.range(0,96).toArray();var model=QPooledPedigreeModelTest.model(d);var geometry=model.geometry(mask);
        var computed=QPooledNullCache.fit(geometry,d.y(),"trait",root,false);assertFalse(computed.reused);assertEquals(1,model.fitExecutions());
        var second=QPooledPedigreeModelTest.model(d);var reused=QPooledNullCache.fit(second.geometry(mask),d.y(),"trait",root,false);assertTrue(reused.reused);assertEquals(0,second.fitExecutions());
        assertArrayEquals(computed.variance,reused.variance,0);assertArrayEquals(computed.py,reused.py,0);
        var rebuilt=QPooledNullCache.fit(geometry,d.y(),"trait",root,true);assertFalse(rebuilt.reused);assertEquals(2,model.fitExecutions());
        assertNotEquals(geometry.nullCacheKey("trait",d.y()),geometry.nullCacheKey("different",d.y()));
        double[] changed=d.y().clone();changed[0]+=.1;assertNotEquals(geometry.nullCacheKey("trait",d.y()),geometry.nullCacheKey("trait",changed));
        Path file;try(var paths=Files.walk(root)){file=paths.filter(p->p.toString().endsWith(".qnull")).findFirst().orElseThrow();}
        byte[] bytes=Files.readAllBytes(file);bytes[100]^=1;Files.write(file,bytes);
        assertThrows(java.io.IOException.class,()->QPooledNullCache.fit(geometry,d.y(),"trait",root,false));
    }
    @Test void batchedNativeRotationsMatchScalarForLargeFamilies(@TempDir Path root)throws Exception {
        Path file=root.resolve("large.tsv");StringBuilder text=new StringBuilder("ID\tSIRE\tDAM\nf0\t0\t0\nm0\t0\t0\nf1\t0\t0\nm1\t0\t0\n");
        int n=80;String[] ids=new String[n],groups=new String[n];double[][] x=new double[n][3];
        for(int i=0;i<n;i++){ids[i]="s"+i;groups[i]=i<40?"FHS":i<60?"WHI":"CARDIA";x[i][0]=1;x[i][1]=i>=40&&i<60?1:0;x[i][2]=i>=60?1:0;if(i<40)text.append(ids[i]+"\tf"+(i/20)+"\tm"+(i/20)+"\n");}
        Files.writeString(file,text);var geometry=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(file),1000).geometry(java.util.stream.IntStream.range(0,n).toArray());
        double[][] values=new double[17][n];for(int a=0;a<17;a++)for(int i=0;i<n;i++)values[a][i]=Math.sin((a+1)*(i+.4));
        double[][] rotated=geometry.rotateRows(values);for(int a=0;a<17;a++)assertArrayEquals(geometry.rotate(values[a]),rotated[a],1e-12);
    }
}
