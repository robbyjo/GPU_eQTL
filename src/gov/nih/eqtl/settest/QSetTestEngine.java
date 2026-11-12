/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import gov.nih.eqtl.*;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.gpu.*;
import gov.nih.jama.QRDecomposition;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.util.*;
import static gov.nih.utils.QDataUtils.kUndefinedValue;
import static gov.nih.eqtl.settest.QSetTestRunner.csv;
import gov.nih.eqtl.settest.QBurdenReference.Variant;
import gov.nih.eqtl.settest.QVariantSetTable.SetDefinition;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;

/** Trait-outer set schedule: one fitted null at a time and atomic result + audit parts. */
public final class QSetTestEngine {
    private QSetTestEngine() { }
    private static final String HEADER = "Set_ID,Trait_ID,Method,Variants,N,DF,Statistic,R_squared,Effect,T,P_value,log10P,P_value_method,Minimum_component_P,Rho_component_P_values,Sample_mask,Null_model";
    private static final String AUDIT = "Set_ID\tRequested\tAbsent\tFrequency_excluded\tIncluded\tStatus\tIncluded_variant_IDs\tDefinition_signature\tMethod\tTrait_ID\tN\tSample_mask\tPattern_monomorphic_or_no_calls\tNull_model";

    public static String validate(QMatrixRowSource traits,int[] columns,double[][] design,QeQTLAnalysisConfig config) throws Exception {
        Set<String> patterns = new HashSet<>();
        MessageDigest content = MessageDigest.getInstance("SHA-256");
        var bits = java.nio.ByteBuffer.allocate(8);
        int skipped = 0, count = 0;
        try (var reader = traits.open(columns)) {
            QMatrixRowSource.Block block;
            while ((block = reader.readBlock(1)) != null) {
                content.update(block.rowIds()[0].getBytes(StandardCharsets.UTF_8)); content.update((byte)0);
                for (double value : block.values()[0]) { bits.clear(); bits.putDouble(value); content.update(bits.array()); }
                int[] observed = QSetTraitModel.observed(block.values()[0]);
                if (observed.length != columns.length && config.getTraitMissingPolicy() != QMissingValuePolicy.PATTERN)
                    throw new IllegalArgumentException("Unfilled trait requires --trait-missing pattern");
                patterns.add(Arrays.toString(observed));
                if (patterns.size() > config.getMaximumTraitPatterns())
                    throw new IllegalArgumentException("Set tests exceed --max-trait-patterns; select an explicit larger limit");
                if (!estimable(design,observed)) {
                    if (!config.getSkipUnestimableTraits()) throw new IllegalArgumentException(
                        "Unestimable set trait mask (rank/DF): "+block.rowIds()[0]+" N="+observed.length);
                    skipped++;
                }
                count++;
            }
        }
        config.getSkatOMethod(); config.getSkatORelativeTolerance(); config.getSkatOMaxEvaluations();
        System.out.println("Set null preflight: traits="+count+"; exact sample masks="+patterns.size()+"; unestimable="+skipped);
        return HexFormat.of().formatHex(content.digest());
    }
    private static boolean estimable(double[][] design,int[] samples) {
        if (samples.length <= design[0].length+1) return false;
        return new QRDecomposition(QSetTraitModel.subsetDesign(design,samples)).isFullRank();
    }

    public static void run(QMatrixRowSource variants,int[] variantColumns,QMatrixRowSource traits,int[] traitColumns,
        double[][] design,QVariantSetTable definitions,QSetTestRunner.Options options,
        QeQTLAnalysisConfig config,QMixedModel mixed,GpuContext[] contexts) throws Exception {
        String traitContent = validate(traits,traitColumns,design,config);
        int setSize = options.setBlockSize() == 0 ? QSetTestRunner.recommendSetBlockSize(definitions,
            variantColumns.length,options.method(),Runtime.getRuntime().maxMemory()/(mixed == null ? 1 : 3)).blockSize() : options.setBlockSize();
        int setBlocks = (definitions.sets().size()+setSize-1)/setSize;
        int traitCount = Math.toIntExact(traits.metadata().rowCount());
        String sourceStamp = QSetTestRunner.signature(variants,variantColumns,traits,traitColumns,design,definitions,options,setSize);
        String modelSources = modelSourceStamp(config);
        String modelSignature = mixed == null ? "fixed" : mixed.signature();
        String extra = config.getTraitMissingPolicy()+";"+config.getFrequencyScope()+";"+config.getSkipUnestimableTraits()
            +";"+config.getSkatOMethod()+";"+config.getSkatORelativeTolerance()+";"+config.getSkatOMaxEvaluations();
        try (GpuContextPool ownership = new GpuContextPool(contexts);
            QSetMatrixOperations operations = new QSetMatrixOperations(contexts,config.get("residualization") != null && config.get("residualization").equals("cpu"))) {
            String signature = hash("set-v4-related-pattern-radial;"+sourceStamp+";"+traitContent+";"+modelSignature+";"+modelSources+";"+extra+";"+operations.executionSignature());
            QAnalysisCheckpoint checkpoint = QAnalysisCheckpoint.open(options.checkpointDirectory(),signature,
                Math.multiplyExact(traitCount,setBlocks),options.resume(),true);
            System.out.println("Set FP64 products: "+operations.executionSignature());
            System.out.println("Backend executes burden aggregation, fixed projection, fitted P*G, genotype scores and score covariance; CPU executes REML, QR, eigenvalues and p-value integration.");
            System.out.println("Set schedule: one trait null at a time; set tile="+setSize+"; related="+(mixed != null)
                +"; frequency scope="+config.getFrequencyScope()+"; SKAT-O="+config.getSkatOMethod());
            try (var reader = traits.open(traitColumns)) {
                for (int trait = 0; trait < traitCount; trait++) {
                    var block = reader.readBlock(1);
                    if (block == null || block.rowOffset() != trait) throw new IOException("Trait source changed");
                    String traitId = block.rowIds()[0];
                    if (traitId.contains("\n") || traitId.contains("\r") || traitId.contains("\t"))
                        throw new IllegalArgumentException("Set trait ID contains a control character");
                    boolean complete = true;
                    for (int tile = 0; tile < setBlocks; tile++) complete &= checkpoint.isComplete(trait*setBlocks+tile);
                    if (complete) continue;
                    int[] samples = QSetTraitModel.observed(block.values()[0]);
                    String mask = hash(Arrays.toString(samples));
                    QSetTraitModel model = estimable(design,samples) ? new QSetTraitModel(block.values()[0],design,mixed) : null;
                    if (model == null && !config.getSkipUnestimableTraits()) throw new IllegalArgumentException("Unestimable set trait "+traitId);
                    for (int tile = 0; tile < setBlocks; tile++) {
                        int part = trait*setBlocks+tile;
                        if (checkpoint.isComplete(part)) continue;
                        QVariantSetTable table = definitions.subset(tile*setSize,Math.min(definitions.sets().size(),(tile+1)*setSize));
                        Map<String,Variant> selected = new LinkedHashMap<>();
                        if (model != null) for (Variant variant : QSetTestRunner.loadSelectedVariants(variants,variantColumns,table,options.traitBlockRows()))
                            selected.put(variant.id(),variant);
                        checkpoint.writeBlock(part,writer -> {
                            for (SetDefinition set : table.sets()) {
                                if (model == null) {
                                    writer.write("A"+audit(set,traitId,samples.length,mask,0,0,List.of(),0,"skipped-unestimable-mask","unfitted",definitions,options));
                                } else analyze(writer,set,selected,model,traitId,mask,definitions,options,config,operations);
                            }
                            if (!modelSources.equals(modelSourceStamp(config)) || !sourceStamp.equals(QSetTestRunner.signature(variants,variantColumns,traits,traitColumns,design,definitions,options,setSize)))
                                throw new IOException("Set input changed before checkpoint commit");
                        });
                    }
                    System.out.println("Set trait "+(trait+1)+"/"+traitCount+": "+traitId+"; N="+samples.length+"; "+(model == null ? "skipped" : model.diagnostics()));
                }
                if (reader.readBlock(1) != null) throw new IOException("Trait source grew during set analysis");
            }
            if (!modelSources.equals(modelSourceStamp(config)) || !sourceStamp.equals(QSetTestRunner.signature(variants,variantColumns,traits,traitColumns,design,definitions,options,setSize)))
                throw new IOException("Set input changed before assembly");
            assemble(checkpoint,options);
            System.out.println("Executed set matrix calls by backend/operation: "+operations.executedProducts());
        }
    }

    private static void analyze(Writer writer,SetDefinition set,Map<String,Variant> variants,QSetTraitModel model,
        String traitId,String mask,QVariantSetTable definitions,QSetTestRunner.Options options,QeQTLAnalysisConfig config,
        QSetMatrixOperations operations) throws Exception {
        int absent = 0, excluded = 0, monomorphic = 0;
        int[] samples = model.samples();
        List<double[]> rows = new ArrayList<>(); List<String> ids = new ArrayList<>();
        QSetTestPolicy policy = options.policy();
        for (var entry : set.entries()) {
            Variant variant = variants.get(entry.variantId());
            if (variant == null) {
                absent++;
                if (policy.absentVariant() == FailurePolicy.ERROR) throw new IllegalArgumentException("Absent set variant "+entry.variantId());
                continue;
            }
            entry.requireExactAlleles(variant.ref(),variant.alt());
            double[] full = variant.altDosagesInternal(), subset = new double[samples.length];
            for (int i = 0; i < samples.length; i++) subset[i] = full[samples[i]];
            double[] frequency = frequency(config.getFrequencyScope().equals("pattern") ? subset : full);
            if (!policy.includes(frequency[0],frequency[1])) { excluded++; continue; }
            double[] observed = frequency(subset);
            if (!Double.isFinite(observed[2])) { monomorphic++; continue; }
            double mean = entry.effectIsAlt() ? observed[2] : 2-observed[2];
            double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < subset.length; i++) {
                if (missing(subset[i])) {
                    if (policy.predictorMissing() == QMissingValuePolicy.ERROR) throw new IllegalArgumentException("Missing set genotype "+variant.id());
                    subset[i] = policy.predictorMissing() == QMissingValuePolicy.MEAN ? mean : 0;
                } else {
                    double dosage = Math.max(0,Math.min(2,subset[i]));
                    subset[i] = entry.effectIsAlt() ? dosage : 2-dosage;
                }
                min = Math.min(min,subset[i]); max = Math.max(max,subset[i]); subset[i] *= entry.weight();
            }
            if (max-min <= 1e-12) { monomorphic++; continue; }
            rows.add(subset); ids.add(variant.id());
        }
        String status = "analyzed", result = null;
        if (rows.isEmpty()) status = "skipped-empty-or-monomorphic";
        else {
            var score = model.scores(rows.toArray(double[][]::new),operations,options.method() == QSetTestMethod.BURDEN);
            double trace = 0; for (int i = 0; i < score.covariance().length; i++) trace += score.covariance()[i][i];
            if (trace < 0 || !Double.isFinite(trace)) throw new ArithmeticException("Invalid set score covariance trace");
            if (!(trace > 0)) status = "skipped-degenerate-kernel";
            else result = result(set.id(),traitId,ids.size(),mask,model,score,options,config);
        }
        if (!status.equals("analyzed") && policy.degenerateSet() == FailurePolicy.ERROR)
            throw new IllegalArgumentException("Set "+set.id()+", trait "+traitId+": "+status);
        writer.write("A"+audit(set,traitId,samples.length,mask,absent,excluded,ids,monomorphic,status,model.diagnostics(),definitions,options));
        if (result != null) writer.write("R"+result);
    }

    private static String result(String set,String trait,int variants,String mask,QSetTraitModel model,
        QSetTraitModel.Scores score,QSetTestRunner.Options options,QeQTLAnalysisConfig config) {
        double statistic,p,logP,rsq = Double.NaN,effect = Double.NaN,t = Double.NaN,minimum = Double.NaN;
        int df = model.nullDf(); String method, components = "";
        if (options.method() == QSetTestMethod.BURDEN) {
            double u = score.u()[0], information = score.covariance()[0][0];
            statistic = u*u/information;
            if (model.related()) {
                df = 1; effect = u/information;
                logP = jdistlib.ChiSquare.cumulative(statistic,1,false,true)/Math.log(10);
                method = "null-reml-score-chi-square";
            } else {
                rsq = statistic/model.nullDf();
                if (rsq > 1+1e-10) throw new ArithmeticException("Invalid burden R-squared");
                rsq = Math.min(1,rsq); df = model.nullDf()-1;
                t = Math.copySign(Math.sqrt(df*rsq/(1-rsq)),u); statistic = t*t;
                effect = u*model.sigma()/information;
                logP = (Math.log(2)+jdistlib.T.cumulative(Math.abs(t),df,false,true))/Math.log(10);
                method = "student-t";
            }
            p = Math.pow(10,logP);
        } else {
            QSetScoreDistribution distribution = new QSetScoreDistribution(score.covariance());
            long seed = options.skatOSeed() ^ set.hashCode() ^ Long.rotateLeft(trait.hashCode(),23);
            var r = options.method() == QSetTestMethod.SKAT ? distribution.skat(score.u())
                : distribution.omnibus(score.u(),options.skatORhoGrid(),config.getSkatOMethod().equals("simulation"),
                    options.skatOSimulations(),seed,config.getSkatORelativeTolerance(),config.getSkatOMaxEvaluations());
            statistic = r.statistic(); p = r.p(); logP = Math.log10(p); method = r.method(); minimum = r.minimumP(); components = r.components();
        }
        boolean keep = switch (options.thresholdType()) {
            case "none" -> true;
            case "pval" -> logP <= Math.log10(options.threshold());
            case "rsq" -> {
                if (model.related() || options.method() != QSetTestMethod.BURDEN) throw new IllegalArgumentException("R-squared threshold requires fixed-effect burden");
                yield rsq >= options.threshold();
            }
            default -> throw new IllegalArgumentException("Invalid set threshold");
        };
        if (!keep) return null;
        return csv(set)+","+csv(trait)+","+options.method().optionName()+","+variants+","+model.sampleCount()+","+df+","+statistic
            +","+rsq+","+effect+","+t+","+p+","+logP+","+csv(method)+","+minimum+","+csv(components)+","+mask+","+csv(model.diagnostics())+System.lineSeparator();
    }
    private static double[] frequency(double[] row) {
        double sum = 0; int count = 0;
        for (double x : row) if (!missing(x)) {
            if (!Double.isFinite(x) || x < -1e-12 || x > 2+1e-12) throw new IllegalArgumentException("Invalid additive set dosage");
            sum += Math.max(0,Math.min(2,x)); count++;
        }
        double mean = sum/count, maf = Math.min(mean/2,1-mean/2);
        return new double[]{maf,2*count*maf,mean};
    }
    private static boolean missing(double value) { return Double.isNaN(value) || value == kUndefinedValue; }
    private static String audit(SetDefinition set,String trait,int n,String mask,int absent,int excluded,List<String> ids,int mono,
        String status,String model,QVariantSetTable definitions,QSetTestRunner.Options options) {
        return set.id()+"\t"+set.entries().size()+"\t"+absent+"\t"+excluded+"\t"+ids.size()+"\t"+status+"\t"+String.join(";",ids)
            +"\t"+definitions.signature()+"\t"+options.method().optionName()+"\t"+trait+"\t"+n+"\t"+mask+"\t"+mono+"\t"+model+System.lineSeparator();
    }
    private static void assemble(QAnalysisCheckpoint checkpoint,QSetTestRunner.Options options) throws Exception {
        Path bundle = checkpoint.directory().resolve("assembled.bundle"); checkpoint.assemble(bundle,"set-bundle-v1");
        Path output = options.output().toAbsolutePath(), audit = options.auditOutput().toAbsolutePath();
        if (output.equals(audit)) throw new IllegalArgumentException("Set output and audit must differ");
        Files.createDirectories(output.getParent()); Files.createDirectories(audit.getParent());
        Path outPart = Files.createTempFile(output.getParent(),"set-result", ".partial");
        Path auditPart = Files.createTempFile(audit.getParent(),"set-audit", ".partial");
        try (var reader = Files.newBufferedReader(bundle); var out = Files.newBufferedWriter(outPart); var aud = Files.newBufferedWriter(auditPart)) {
            reader.readLine(); out.write(HEADER); out.newLine(); aud.write(AUDIT); aud.newLine();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("R")) { out.write(line.substring(1)); out.newLine(); }
                else if (line.startsWith("A")) { aud.write(line.substring(1)); aud.newLine(); }
                else throw new IOException("Invalid set checkpoint record");
            }
        }
        move(auditPart,audit); move(outPart,output); Files.delete(bundle);
        if (!options.keepCheckpoints()) checkpoint.removeCompleted();
    }
    private static void move(Path from,Path to) throws IOException {
        try { Files.move(from,to,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING); }
        catch (AtomicMoveNotSupportedException e) { Files.move(from,to,StandardCopyOption.REPLACE_EXISTING); }
    }
    private static String hash(String value) throws Exception {
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value.getBytes(StandardCharsets.UTF_8)));
    }
    private static String modelSourceStamp(QeQTLAnalysisConfig config) throws IOException {
        StringBuilder result = new StringBuilder();
        for (String file : new String[]{config.getCovariateFilename(),config.getPedigreeFilename(),config.getRelationshipMatrixFilename()}) {
            if (file == null) continue;
            Path path = Path.of(file).toAbsolutePath().normalize();
            var a = Files.readAttributes(path,java.nio.file.attribute.BasicFileAttributes.class);
            result.append(path).append(':').append(a.size()).append(':').append(a.lastModifiedTime()).append(':').append(a.fileKey()).append('\n');
        }
        return result.toString();
    }
}
