/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.List;

/** Thread-safe pool that gives each worker exclusive access to a GPU context. */
public final class GpuContextPool implements AutoCloseable {
	private final List<GpuContext> allContexts;
	private final Deque<GpuContext> availableContexts;
	private boolean closed;

	public GpuContextPool(GpuContext[] contexts) {
		this(contexts == null ? Collections.<GpuContext>emptyList() : Arrays.asList(contexts));
	}

	public GpuContextPool(Collection<? extends GpuContext> contexts) {
		if (contexts == null) {
			throw new IllegalArgumentException("contexts must not be null");
		}
		allContexts = new ArrayList<GpuContext>(contexts);
		availableContexts = new ArrayDeque<GpuContext>(contexts);
	}

	public synchronized List<GpuContext> getAllContexts() {
		return Collections.unmodifiableList(new ArrayList<GpuContext>(allContexts));
	}

	public synchronized int getNumContexts() {
		return allContexts.size();
	}

	public synchronized GpuContext reserveContext() {
		while (availableContexts.isEmpty() && !closed) {
			try {
				wait();
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new GpuException("Interrupted while waiting for a GPU context", e);
			}
		}
		if (closed) {
			throw new GpuException("GPU context pool is closed");
		}
		return availableContexts.removeFirst();
	}

	public synchronized void releaseContext(GpuContext context) {
		if (!closed && allContexts.contains(context) && !availableContexts.contains(context)) {
			availableContexts.addLast(context);
			notifyAll();
		}
	}

	@Override
	public synchronized void close() {
		if (closed) {
			return;
		}
		closed = true;
		for (GpuContext context : allContexts) {
			context.close();
		}
		availableContexts.clear();
		notifyAll();
	}
}
