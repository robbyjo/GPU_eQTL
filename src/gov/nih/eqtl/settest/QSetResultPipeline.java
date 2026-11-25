/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.io.Writer;
import java.util.ArrayDeque;
import java.util.concurrent.*;

/** Bounded CPU-only tails, consumed in source order by the context owner. */
final class QSetResultPipeline implements AutoCloseable {
    private final ExecutorService executor;
    private final ArrayDeque<Future<String>> pending = new ArrayDeque<>();
    private final int workers;

    static int workers(int requested,int sets,int cores,long heap) {
        if (requested < 0 || sets < 1 || cores < 1 || heap < 1)
            throw new IllegalArgumentException("Invalid set worker budget");
        int desired = requested == 0 ? Math.min(8,cores) : requested;
        return (int)Math.max(1,Math.min(Math.min(desired,Math.min(sets,cores)),
            Math.min(256,heap/(32L*1024*1024))));
    }
    QSetResultPipeline(int workers) {
        if (workers < 1) throw new IllegalArgumentException("Positive set worker count required");
        this.workers=workers;
        executor=Executors.newFixedThreadPool(workers,r -> {
            Thread thread=new Thread(r,"set-tail-worker");thread.setDaemon(true);return thread;
        });
    }
    void submit(Writer writer,Callable<String> task) throws Exception {
        // At most workers compact states are retained. Never queue G, P*G, C,
        // native algebra or a GPU context; those stay on the producer thread.
        if (pending.size() == workers) writeFirst(writer);
        pending.addLast(executor.submit(task));
    }
    void drain(Writer writer) throws Exception {
        while (!pending.isEmpty()) writeFirst(writer);
    }
    private void writeFirst(Writer writer) throws Exception {
        try { writer.write(pending.getFirst().get());pending.removeFirst(); }
        catch (ExecutionException failure) {
            Throwable cause=failure.getCause();
            if (cause instanceof Exception exception) throw exception;
            if (cause instanceof Error error) throw error;
            throw new IllegalStateException(cause);
        }
    }
    @Override public void close() {
        for (Future<String> task:pending) task.cancel(true);
        executor.shutdownNow();
        boolean interrupted=false;
        for (;;) try {
            if (executor.awaitTermination(1,TimeUnit.SECONDS)) break;
        } catch (InterruptedException failure) { interrupted=true; }
        if (interrupted) Thread.currentThread().interrupt();
    }
}
