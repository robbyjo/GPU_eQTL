/*
 * Roby Joehanes
 * 
 * Copyright 2009 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package gov.nih.eqtl.datastructure;

import static gov.nih.eqtl.datastructure.EBase.*;
import gov.nih.utils.QDataUtils;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QSNPDataInt extends QSNPData {
	/**
	 * SNP data are stored into n x 2 bits storage. Each SNP is represented by: 1, 2, 3, or 0. 1 and 3 are homozygotes,
	 * while 2 is heterozygote. 0 means missing.
	 * 
	 */
	protected long[] mSNPData;
	protected int
		mMajorAlleleCount,
		mMinorAlleleCount;
	protected EBase mMajorAllele, mMinorAllele;
	protected byte[] mRawGenotypes;

	/*
	protected int
		mMajorHomozygote,
		mMinorHomozygote,
		mHeterozygote,
		mMajorHomozygoteCount,
		mMinorHomozygoteCount,
		mHeterozygoteCount;
	protected int[] mAlleleCounts;
	protected int[] mDiploidAlleleCounts;
	//*/

	static final boolean recordRawGenotypes = false;

	// == for (int i = 0; i < 32; i++) mMask[31-i] = ~(3L << (i << 1));
	// for (int i = 0; i < mMask.length; i++)
	//    System.out.print("0x"+Long.toHexString(mMask[i]).toUpperCase() + "L, ");
	protected static final long[] mMask = new long[] {
		0x3FFFFFFFFFFFFFFFL, 0xCFFFFFFFFFFFFFFFL, 0xF3FFFFFFFFFFFFFFL, 0xFCFFFFFFFFFFFFFFL,
		0xFF3FFFFFFFFFFFFFL, 0xFFCFFFFFFFFFFFFFL, 0xFFF3FFFFFFFFFFFFL, 0xFFFCFFFFFFFFFFFFL,
		0xFFFF3FFFFFFFFFFFL, 0xFFFFCFFFFFFFFFFFL, 0xFFFFF3FFFFFFFFFFL, 0xFFFFFCFFFFFFFFFFL,
		0xFFFFFF3FFFFFFFFFL, 0xFFFFFFCFFFFFFFFFL, 0xFFFFFFF3FFFFFFFFL, 0xFFFFFFFCFFFFFFFFL,
		0xFFFFFFFF3FFFFFFFL, 0xFFFFFFFFCFFFFFFFL, 0xFFFFFFFFF3FFFFFFL, 0xFFFFFFFFFCFFFFFFL,
		0xFFFFFFFFFF3FFFFFL, 0xFFFFFFFFFFCFFFFFL, 0xFFFFFFFFFFF3FFFFL, 0xFFFFFFFFFFFCFFFFL,
		0xFFFFFFFFFFFF3FFFL, 0xFFFFFFFFFFFFCFFFL, 0xFFFFFFFFFFFFF3FFL, 0xFFFFFFFFFFFFFCFFL,
		0xFFFFFFFFFFFFFF3FL, 0xFFFFFFFFFFFFFFCFL, 0xFFFFFFFFFFFFFFF3L, 0xFFFFFFFFFFFFFFFCL
	};
	protected static final long[] mLeftMask = new long[] {
		0x0000000000000000L, 0xC000000000000000L, 0xF000000000000000L, 0xFC00000000000000L,
		0xFF00000000000000L, 0xFFC0000000000000L, 0xFFF0000000000000L, 0xFFFC000000000000L,
		0xFFFF000000000000L, 0xFFFFC00000000000L, 0xFFFFF00000000000L, 0xFFFFFC0000000000L,
		0xFFFFFF0000000000L, 0xFFFFFFC000000000L, 0xFFFFFFF000000000L, 0xFFFFFFFC00000000L,
		0xFFFFFFFF00000000L, 0xFFFFFFFFC0000000L, 0xFFFFFFFFF0000000L, 0xFFFFFFFFFC000000L,
		0xFFFFFFFFFF000000L, 0xFFFFFFFFFFC00000L, 0xFFFFFFFFFFF00000L, 0xFFFFFFFFFFFC0000L,
		0xFFFFFFFFFFFF0000L, 0xFFFFFFFFFFFFC000L, 0xFFFFFFFFFFFFF000L, 0xFFFFFFFFFFFFFC00L,
		0xFFFFFFFFFFFFFF00L, 0xFFFFFFFFFFFFFFC0L, 0xFFFFFFFFFFFFFFF0L, 0xFFFFFFFFFFFFFFFCL
	};
	protected static final long[] mRightMask = new long[] {
		0x3FFFFFFFFFFFFFFFL, 0x0FFFFFFFFFFFFFFFL, 0x03FFFFFFFFFFFFFFL, 0x00FFFFFFFFFFFFFFL,
		0x003FFFFFFFFFFFFFL, 0x000FFFFFFFFFFFFFL, 0x0003FFFFFFFFFFFFL, 0x0000FFFFFFFFFFFFL,
		0x00003FFFFFFFFFFFL, 0x00000FFFFFFFFFFFL, 0x000003FFFFFFFFFFL, 0x000000FFFFFFFFFFL,
		0x0000003FFFFFFFFFL, 0x0000000FFFFFFFFFL, 0x00000003FFFFFFFFL, 0x00000000FFFFFFFFL,
		0x000000003FFFFFFFL, 0x000000000FFFFFFFL, 0x0000000003FFFFFFL, 0x0000000000FFFFFFL,
		0x00000000003FFFFFL, 0x00000000000FFFFFL, 0x000000000003FFFFL, 0x000000000000FFFFL,
		0x0000000000003FFFL, 0x0000000000000FFFL, 0x00000000000003FFL, 0x00000000000000FFL,
		0x000000000000003FL, 0x000000000000000FL, 0x0000000000000003L, 0x0000000000000000L
	};
	protected static final int[][] mDiploidIdx = new int[][] {
		{0, 1, 2, 3}, {-1, 4, 5, 6}, {-1, -1, 7, 8}, {-1, -1, -1, 9}
	};
	protected static final String[]
		sAlleleSymbols = new String[] { A.name(), C.name(), G.name(), T.name()},
		sDiploidSymbols = new String[] { A.pair(A), A.pair(C), A.pair(G), A.pair(T),
			C.pair(C), C.pair(G), C.pair(T),
			G.pair(G), G.pair(T), T.pair(T)};

	public QSNPDataInt() {}

	public QSNPDataInt(String id, String[] tokens)
	{
		setID(id);
		setSNPValues(tokens);
	}

	private static final int decodeAllele(String token)
	{
		char cc = token.charAt(0);
		switch (cc) {
			case 'A': return A.ordinal();
			case 'C': return C.ordinal();
			case 'G': return G.ordinal();
			case 'T': return T.ordinal();
			case '-':
			case '?':
			case '0': return kMissingCode;
		}
		System.err.println(String.format("Unexpected token '%s'. Considered as missing.", token));
		return kMissingCode;
	}

	@Override
	public void setSNPValues(String[] tokens)
	{
		int numTokens = tokens.length;
		mNumIndividuals = numTokens / 2;

		byte[] genotypes = new byte[mNumIndividuals];
		if (recordRawGenotypes)
			mRawGenotypes = new byte[mNumIndividuals];
		int[]
			genotypeCount = new int[10],
			alleleCount = new int[4],
			allelePos = new int[] {0,1,2,3};
		//mDiploidAlleleCounts = genotypeCount;
		for (int i = 0; i < numTokens; i++)
		{
			int
				t1 = decodeAllele(tokens[i++]),
				t2 = decodeAllele(tokens[i]);
			if (recordRawGenotypes)
				mRawGenotypes[i/2] = (byte) ((t1+1) * 5 + (t2+1));
			if (t2 < t1)
			{
				int temp = t2; t2 = t1; t1 = temp;
			}

			//-t1*t1/2 + 4.5*t1 + t2 - t1
			if (t1 != kMissingCode)
			{
				alleleCount[t1]++;
				alleleCount[t2]++;
				int genotype = genotypes[i>>>1] = (byte) mDiploidIdx[t1][t2];
				genotypeCount[genotype]++;
			} else genotypes[i>>>1] = kMissingCode;
		}
		//mAlleleCounts = new int[4]; System.arraycopy(alleleCount, 0, mAlleleCounts, 0, 4);
		QDataUtils.sort2DByRow(new int[][] {alleleCount, allelePos}, 0);
		int
			majorAllele = allelePos[3],
			minorAllele = allelePos[2];
		if (alleleCount[1] > 0 || alleleCount[0] > 0)
		{
			System.err.print(String.format("Major allele %s, minor allele %s. ", sAlleleSymbols[majorAllele], sAlleleSymbols[minorAllele]));
			System.err.println(String.format("However, allele %s is present with %d (allele %s = %d), which will be considered as missing.",
				sAlleleSymbols[allelePos[1]], alleleCount[1], sAlleleSymbols[allelePos[0]], alleleCount[0]));
		}

		mMajorAllele = EBase.values()[majorAllele];
		mMinorAllele = EBase.values()[minorAllele];
		mMajorAlleleCount = alleleCount[3];
		mMinorAlleleCount = alleleCount[2];
		int majorHomozygote = mDiploidIdx[majorAllele][majorAllele];
		int minorHomozygote = mDiploidIdx[minorAllele][minorAllele];
		int heterozygote = majorAllele < minorAllele ? mDiploidIdx[majorAllele][minorAllele] : mDiploidIdx[minorAllele][majorAllele];

		mSNPData = new long[mNumIndividuals/32 + ((mNumIndividuals & 31) > 0 ? 1 : 0)];
		for (int i = 0; i < mNumIndividuals; i++)
		{
			int genotype = genotypes[i];
			int code = genotype == majorHomozygote ? 0 : genotype == heterozygote ? 1 : genotype == minorHomozygote ? 2 : kMissingCode;
			if (code != kMissingCode) {
				//mSNPData[longIdx] = (mSNPData[longIdx] & mMask[bitOffset]) | (((long) (code+1)) << ((31-bitOffset) << 1));
				setSNPCode(mSNPData, i, code);
				//System.out.println(Long.toHexString(mSNPData[longIdx]));
				/*
				switch (code) {
					case 0: mMajorHomozygoteCount++; break;
					case 1: mHeterozygoteCount++; break;
					case 2: mMinorHomozygoteCount++;
				}
				//*/
			} //else mNumMissing++;
		}
		/*
		for (int i = 0; i < mNumIndividuals; i++) {
			String t1 = tokens[2*i];
			String t2 = tokens[2*i+1];
			System.out.println(t1 + t2 + "=" + Long.toHexString(mSNPData[i >>> 5]) + "=" + getSNPCode(i));
		}
		//*/
	}

	/**
	 * Encoded version
	 * @param tokens
	 */
	public void setSNPValues(int[] tokens)
	{
		mNumIndividuals = tokens.length;
		mSNPData = new long[mNumIndividuals/32 + ((mNumIndividuals & 31) > 0 ? 1 : 0)];
		for (int i = 0; i < mNumIndividuals; i++)
		{
			int code = tokens[i];
			if (code != kMissingCode)
				setSNPCode(mSNPData, i, code);
		}
	}

	private static final int getSNPCode(long[] snpData, int subjectNo)
	{
		int offset = subjectNo & 31;
		return (int) ((snpData[subjectNo >>> 5] & ~mMask[offset]) >>> ((31 - offset) << 1)) - 1;
	}

	private static final void setSNPCode(long[] snpData, int subjectNo, int code)
	{
		code = (code+1) & 3;
		int
			bitOffset = subjectNo & 31,
			longIdx = subjectNo >>> 5;
		snpData[longIdx] = (snpData[longIdx] & mMask[bitOffset]) | (((long) code) << ((31-bitOffset) << 1));
	}

	public int getSNPCode(int subjectNo)
	{	return getSNPCode(mSNPData, subjectNo); }

	/**
	 * Set SNP code for a given subject. Code MUST be one of:<br>
	 * 0 = major homozygote<br>
	 * 1 = heterozygote<br>
	 * 2 = minor homozygote
	 * -1 = missing<br>
	 * @param subjectNo
	 * @param code
	 */
	public void setSNPCode(int subjectNo, int code)
	{	setSNPCode(mSNPData, subjectNo, code); }

	public void removeSubject(int subjectNo)
	{
		int
			bitOffset = subjectNo & 31,
			longIdx = subjectNo >>> 5,
			snpArraylengthMin1 = mSNPData.length - 1;
		long nextCode = longIdx < snpArraylengthMin1? (mSNPData[longIdx+1] & ~mMask[0]) >>> 62 : 0;
		mSNPData[longIdx] = (mSNPData[longIdx] & mLeftMask[bitOffset]) | ((mSNPData[longIdx] & mRightMask[bitOffset]) << 2) | nextCode;
		for (int i = longIdx + 1; i < snpArraylengthMin1; i++)
			mSNPData[i] = (mSNPData[i] << 2) | ((mSNPData[i+1] & ~mMask[0]) >>> 62);
		if (longIdx < snpArraylengthMin1)
			mSNPData[snpArraylengthMin1] <<= 2;
		if (recordRawGenotypes) {
			byte[] newRawGenotypes = new byte[mNumIndividuals - 1];
			System.arraycopy(mRawGenotypes, 0, newRawGenotypes, 0, subjectNo);
			System.arraycopy(mRawGenotypes, subjectNo+1, newRawGenotypes, subjectNo, mNumIndividuals - subjectNo - 1);
			mRawGenotypes = newRawGenotypes;
		}
		mNumIndividuals--;
	}

	public void rearrangeSubjects(int[] newIdx)
	{
		int
			newNumInds = 0,
			j = 0,
			indexLen = newIdx.length;
		for (int i = 0; i < indexLen; i++)
			if (newIdx[i] >= 0) newNumInds++;
		long[] newSNPData = new long[newNumInds/32 + ((newNumInds & 31) > 0 ? 1 : 0)];
		byte[] newRawGenotypes = null;
		if (recordRawGenotypes)
			newRawGenotypes = new byte[newNumInds];
		for (int i = 0; i < indexLen; i++)
		{
			if (newIdx[i] >= 0)
			{
				setSNPCode(newSNPData, j, getSNPCode(mSNPData, newIdx[i]));
				if (recordRawGenotypes)
					newRawGenotypes[j] = mRawGenotypes[newIdx[i]];
				j++;
			}
		}
		if (recordRawGenotypes)
			mRawGenotypes = newRawGenotypes;
		mSNPData = newSNPData;
		mNumIndividuals = newNumInds;
	}

	@Override
	public int[] getSNPCodes()
	{
		int[] codes = new int[mNumIndividuals];
		for (int i = 0; i < mNumIndividuals; i++)
			codes[i] = getSNPCode(i);
		return codes;
	}

	/**
	 * Same as getSNPCodes, but returns double values
	 * @return
	 */
	@Override
	public double[] getSNPValues()
	{
		double[] codes = new double[mNumIndividuals];
		for (int i = 0; i < mNumIndividuals; i++)
			codes[i] = getSNPCode(i);
		return codes;
	}

	public String[] getRawGenotypes()
	{
		String[] rawGenotypes = new String[mNumIndividuals];
		String[] codes = new String[] { "0", A.name(), C.name(), G.name(), T.name() };
		for (int i = 0; i < mNumIndividuals; i++) {
			int
				t1 = mRawGenotypes[i],
				t2 = t1 % 5;
			t1 = t1 / 5;
			rawGenotypes[i] = codes[t1] + codes[t2];
		}
		return rawGenotypes;
	}

	/*
	public int[] getAlleleCounts()
	{	return mAlleleCounts; }

	public int[] getDiploidAlleleCounts()
	{	return mDiploidAlleleCounts; }

	public int getNumDiploidAlleleTypes()
	{
		int ct = 0;
		for (int count: mDiploidAlleleCounts)
			if (count > 0) ct++;
		return ct;
	}

	public String getMajorHomozygoteAllele()
	{	return sDiploidSymbols[mMajorHomozygote]; }

	public String getMinorHomozygoteAllele()
	{	return sDiploidSymbols[mMinorHomozygote]; }

	public String getHeterozygoteAllele()
	{	return sDiploidSymbols[mHeterozygote]; }

	public int getMajorHomozygoteCount()
	{	return mMajorHomozygoteCount; }

	public int getMinorHomozygoteCount()
	{	return mMinorHomozygoteCount; }

	public int getHeterozygoteCount()
	{	return mHeterozygoteCount; }
	//*/

	@Override
	public int getNumMissing()
	{
		int numMissing = 0;
		for (int i = 0; i < mNumIndividuals; i++)
			if (getSNPCode(mSNPData, i) == kMissingCode)
				numMissing++;

		return numMissing;
	}

	public double getMinorAlleleFreq()
	{	return mMinorAlleleCount * 1.0 / (mMinorAlleleCount + mMajorAlleleCount); }

	@Override
	public void purge()
	{
		mSNPData = null;
		mRawGenotypes = null;
	}
}
