/*
 * Roby Joehanes
 * 
 * Copyright 2007 Roby Joehanes
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
/**
 * 
 */
package gov.nih.utils;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.util.StringTokenizer;

/**
 * Simple, hand-made lexer class
 * 
 * @author Roby Joehanes
 *
 */
public class QDefaultLexer
{
	public static final String
		sComment1 = "#", //$NON-NLS-1$
		sComment2 = "//", //$NON-NLS-1$
		sUTF8BOM = "﻿", //$NON-NLS-1$
		sUTF16BOMBigEndian = "\uFFFE", //$NON-NLS-1$
		sUTF16BOMLittleEndian = "\uFEFF"; //$NON-NLS-1$

	protected LineNumberReader mFileReader;
	protected String mPeekedToken;
	protected String mCurrentLine;
	protected StringTokenizer mTokenizer;
	protected int mColumn = 0;
	protected int mLastTokenLength = 0;

	public QDefaultLexer(Reader reader)
	{	mFileReader = new LineNumberReader(reader);	}

	public int getLineNumber()
	{	return mFileReader.getLineNumber(); }

	public int getColumn()
	{	return mColumn + 1; }

	public boolean matchNextToken(String str) throws IOException
	{
		String token = readNextToken();
		if (token == null)
			return false;
		return token.equalsIgnoreCase(str);
	}

	public boolean matchPeekedToken(String str) throws IOException
	{
		String token = peekToken();
		if (token == null)
			return false;
		return token.equalsIgnoreCase(str);
	}

	public String readNextToken() throws IOException
	{
		peekToken();
		String peekedToken = mPeekedToken;
		mColumn += mLastTokenLength;
		if (peekedToken != null)
		{
			mLastTokenLength = peekedToken.length();
			mColumn = mCurrentLine.indexOf(peekedToken, mColumn);
		}
		else
			mLastTokenLength = 0;
		mPeekedToken = null;
		return peekedToken;
	}

	public boolean reachedEOL()
	{	return !mTokenizer.hasMoreTokens(); }

	public String peekToken() throws IOException
	{
		// If we already peeked the token, just return the token,
		// otherwise read it.
		if (mPeekedToken != null)
			return mPeekedToken;
		if (mTokenizer == null || !mTokenizer.hasMoreTokens())
		{
			mTokenizer = null; // Tell Java to garbage-collect this tokenizer
			String line = readLine();
			if (line == null)
				return null;
			mTokenizer = new StringTokenizer(line);
		}
		mPeekedToken = mTokenizer.nextToken();
		if (mPeekedToken.equals(sUTF8BOM) || mPeekedToken.equals(sUTF16BOMLittleEndian) || mPeekedToken.equals(sUTF16BOMBigEndian))
		{
			mPeekedToken = null;
			return peekToken();
		}
		return mPeekedToken;
	}

	public String readUntilEOL()
	{
		if (mTokenizer == null)
			return null;
		// mTokenizer must be non-null
		StringBuffer buffer = new StringBuffer();
		if (mPeekedToken != null)
		{
			buffer.append(mPeekedToken);
			mPeekedToken = null;
		}
		while (mTokenizer.hasMoreTokens())
			buffer.append(" " + mTokenizer.nextToken()); //$NON-NLS-1$
		return buffer.toString().trim();
	}

	public boolean isEOF() throws IOException
	{	return peekToken() == null; }

	public String readLine() throws IOException
	{
		do
		{
			mCurrentLine = mFileReader.readLine();
			mColumn = 0;
			mLastTokenLength = 0;
			if (mCurrentLine != null)
			{
				mCurrentLine = mCurrentLine.trim();
				int commentIndex = mCurrentLine.indexOf(sComment1);
				int commentIndex2 = mCurrentLine.indexOf(sComment2);
				if (commentIndex != -1)
				{
					if (commentIndex2 > -1 && commentIndex2 < commentIndex)
						commentIndex = commentIndex2;
				}
				else
					commentIndex = commentIndex2;
				if (commentIndex > -1) // Chop the comment off (if any)
					mCurrentLine = mCurrentLine.substring(0, commentIndex).trim();
				if (mCurrentLine.length() == 0)
					continue;
				return mCurrentLine;
			}
			return null;
		} while (true);
	}
}