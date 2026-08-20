/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
package gov.nih.eqtl;

import static java.lang.Math.log;
import static java.lang.Math.log10;
import static java.lang.Math.sqrt;

import net.sourceforge.jdistlib.T;

/** Production eQTL statistics derived from the standardized matrix product. */
public final class QeQTLStatistics {
    private static final double LOG10_OF_2 = log10(2);
    private static final double LOG_OF_10 = log(10);

    public record Result(double rSquared, double effect, double tStatistic, double log10P) { }

    private QeQTLStatistics() { }

    public static Result calculate(double correlation, double expressionStandardDeviation,
        double genotypeStandardDeviation, int errorDegreesOfFreedom, int degreesOfFreedomOffset) {
        double rSquared = correlation * correlation;
        double effect = correlation * expressionStandardDeviation / genotypeStandardDeviation;
        double absoluteT = sqrt(rSquared * errorDegreesOfFreedom / (1 - rSquared));
        double log10P = LOG10_OF_2 + T.cumulative(absoluteT,
            errorDegreesOfFreedom - degreesOfFreedomOffset, false, true) / LOG_OF_10;
        return new Result(rSquared, effect, Math.copySign(absoluteT, correlation), log10P);
    }
}
