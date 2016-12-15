/*
 * Roby Joehanes
 * 
 * Copyright 2013 Roby Joehanes
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
package gov.nih.tools;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.tools.FileObject;
import javax.tools.ForwardingJavaFileManager;
import javax.tools.JavaCompiler;
import javax.tools.JavaFileManager;
import javax.tools.JavaFileObject.Kind;
import javax.tools.SimpleJavaFileObject;
import javax.tools.ToolProvider;

import com.csvreader.CsvReader;

import static java.lang.Math.*;

/**
 * Very powerful file lookup tool
 * @author Roby Joehanes
 */
public class FilterLines {

	// java -Xmx64G -cp .:javacsv.jar FilterLines "Rs_ID %% file('bpsnps-20130725.csv',8) || replaceAll(Marker,':[ACGT].*', '') %% file('bpsnps-20130725.csv',7)" eqtl-gene-1000g-result-p1e-4-reannot-irsq1-maf01-fdr.txt bpsnps-gene.txt
	public static final void main(String[] args) {
		if (ToolProvider.getSystemJavaCompiler() == null) {
			System.err.println("Error! Compiler subsystem is not installed correctly! ToolProvider.getSystemJavaCompiler() returns null! Program aborted!");
			System.exit(-1);
		}
		boolean countOnly = args[0].equals("--count");
		if (countOnly) {
			args[0] = args[1];
			args[1] = args[2];
		}

		System.out.println("Expression: " + args[0]);
		TokenManager tokmgr = TokenManager.tokenize(args[0]);
		ASTNode node = ASTBuilder.parseExpr(tokmgr);
		if (node == null)
			throw new RuntimeException("Parse error!");
		ASTBuilder.inferTypes(node);
		List<ASTNode> termNodes = ASTBuilder.getTerminalNodes(node);
		Method method = null;

		String
			inputfile = args[1],
			outputfile = countOnly ? null : args[2];
		System.out.println("Input file: " + inputfile);
		System.out.println("Output file: " + outputfile);
		CsvReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0, lineCut = 0;
		try {
			reader = new CsvReader(inputfile);
			if (!countOnly)
				writer = new PrintWriter(outputfile);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					if (writer != null)
						writer.println(reader.getRawRecord());
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i], i);
					IdentityHashMap<ASTNode, String> symbolTable = new IdentityHashMap<ASTNode, String>();
					for (ASTNode curNode: termNodes)
						if (curNode.getNodeType() == ASTNodeType.IDENTIFIER) {
							String tok = curNode.getToken();
							Integer idx = colnamesToIdx.get(tok);
							if (idx == null) {
								if (writer != null) writer.close(); writer = null;
								reader.close(); reader = null;
								throw new RuntimeException("Column " + tok + " is not found!");
							}
							symbolTable.put(curNode, idx.toString());
						}
					// Generate a Java code for a filter
					String code = CodeGenerator.generate(node, symbolTable);
					System.out.println("Generated code:");
					System.out.println(code);
					// Compile and load on the fly
					MemoryClassLoader mcl = new MemoryClassLoader(CodeGenerator.className, code);
					method = mcl.loadClass(CodeGenerator.className).getMethod(CodeGenerator.methodName, String[].class);
					continue;
				}
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
				try {
					if (((Boolean) method.invoke(null, (Object) tokens)).booleanValue()) {
						if (writer != null) {
							writer.print(tokens[0]);
							for (int i = 1; i < tokens.length; i++)
								writer.print("," + tokens[i]);
							writer.println();
						}
						lineCut++;
					}
				} catch (Exception exx) {
					// Most likely due to missing value
				}
			};
			System.out.println(lineNo + " lines were read.");
			System.out.println(lineCut + " lines were written.");
			if (writer != null) writer.close(); writer = null;
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (Exception e) {}
			if (writer != null)
				try {
					writer.close();
				} catch (Exception e) {}
		}
	}
}

class TokenManager {
	public int[] token_types;
	public String[] tokens;
	private int pos = 0;

	public TokenManager(String[] t, int[] types) {
		tokens = t; token_types = types;
	}

	// This is hand-made lexer. I did only minimum error check!
	static final int
		state_invalid = -1,
		state_clear = 0,
		state_num = 1,
		state_long = 2,
		state_float = 3,
		state_float_e = 4,
		state_float_e_noplusmin = 5,
		state_float_e_plusmin = 6,
		state_ident = 7,
		state_oper = 8,
		state_oper_dbl = 9,
		state_oper_tpl = 10,
		state_oper_compeql = 11,
		state_string = 12,
		state_string2 = 13,
		state_lookup[] = new int[] { state_clear, state_long, state_long, state_float, state_float, state_float, state_float,
			state_ident, state_oper, state_oper, state_oper, state_oper, state_string, state_string };

	static Set<String> acceptedOperatorSet = new HashSet<String>();
	static Set<String> numericFunCall = new HashSet<String>();
	static Set<String> stringFunCall = new HashSet<String>();
	static Set<String> stringReturningIntFunCall = new HashSet<String>();
	static Set<String> stringReturningBoolFunCall = new HashSet<String>();
	static {
		// %% means "in"
		final String[] operSet = new String[] { "+", "-", "*", "**", "/", "%", "%%", "(", ")", "[", "]", "{", "}",
			"<", "<=", ">", ">=", "==", "<<", ">>", ">>>", "&", "&&", "|", "||", "^", "^^", "~", "!", ",", "?", ":"};
		for (String oper: operSet)
			acceptedOperatorSet.add(oper);
		String[] funSet = new String[] { "abs", "log", "log10", "exp", "pow", "round", "floor", "ceil", "atan", "asin", "acos",
			"atan2", "cbrt", "sqrt", "cos", "sin", "tan", "anh", "cosh", "sinh", "expm1", "log1p", "Double.parseDouble",
			"Integer.parseInt", "Long.parseLong"};
		for (String fun: funSet)
			numericFunCall.add(fun);
		funSet = new String[] { "upper", "lower", "toUpperCase", "toLowerCase", "trim", "substring", "replaceAll", "replaceFirst" };
		for (String fun: funSet)
			stringFunCall.add(fun);
		funSet = new String[] { "indexOf", "lastIndexOf" };
		for (String fun: funSet)
			stringReturningIntFunCall.add(fun);
		funSet = new String[] { "startsWith", "endsWith", "contains", "matches" };
		for (String fun: funSet)
			stringReturningBoolFunCall.add(fun);
	}

	public static final TokenManager tokenize(String line) {
		List<String> tokens = new ArrayList<String>();
		char[] chars = line.toCharArray();
		int len = chars.length;
		StringBuilder buf = new StringBuilder();
		List<Integer> states = new ArrayList<Integer>();
		int state = state_clear;
		char lastOper = 0;
		for (int i = 0; i < len; i++) {
			char ch = chars[i];

			// Space means clearing token
			if (ch == ' ') {
				if (buf.length() > 0) {
					tokens.add(buf.toString());
					buf.delete(0, buf.length());
					states.add(state_lookup[state]);
				}
				state = state_clear;
				continue;
			}
			boolean accepted = false;
			switch(state) {
				case state_num:
					if ('0' <= ch && ch <= '9') {
						buf.append(ch);
						accepted = true;
					} else if (ch == 'l' || ch == 'L') {
						accepted = true;
						continue;
					} else if (ch == '.') {
						buf.append(ch);
						state = state_float;
						accepted = true;
					} else if (ch == 'e' || ch == 'E') {
						buf.append(ch);
						state = state_float_e;
						accepted = true;
					}
					break;
				case state_float:
					if ('0' <= ch && ch <= '9') {
						buf.append(ch);
						accepted = true;
					} else if (ch == 'l' || ch == 'L' || ch == '.') {
						throw new RuntimeException();
					} else if (ch == 'e' || ch == 'E') {
						buf.append(ch);
						state = state_float_e;
						accepted = true;
					}
					break;
				case state_float_e:
					if ('0' <= ch && ch <= '9') {
						buf.append(ch);
						state = state_float_e_noplusmin;
						accepted = true;
					} else if (ch == '+' || ch == '-') {
						buf.append(ch);
						state = state_float_e_plusmin;
						accepted = true;
					}
					break;
				case state_float_e_noplusmin:
					if ('0' <= ch && ch <= '9') {
						buf.append(ch);
						accepted = true;
					}
					break;
				case state_float_e_plusmin:
					if ('0' <= ch && ch <= '9') {
						buf.append(ch);
						accepted = true;
					} else if (ch == '+' || ch == '-') {
						state = state_clear;
						accepted = false;
					}
					break;
				case state_ident:
					if (('a' <= ch && ch <= 'z') || ('A' <= ch && ch <= 'Z') || ('0' <= ch && ch <= '9') || ch == '.' || ch == '_') {
						buf.append(ch);
						accepted = true;
					}
					break;
				case state_oper:
					switch (ch) {
						case '+':
						case '-':
						case '/':
						case '~':
						case '!':
						case ',':
						case '(':
						case ')':
						case '[':
						case ']':
						case '{':
						case '}':
						case '?':
						case ':':
							buf.append(ch);
							accepted = true;
							tokens.add(buf.toString());
							buf.delete(0, buf.length());
							states.add(state_lookup[state]);
							state = state_clear;
							lastOper = 0;
							break;
						case '*':
						case '%':
						case '|':
						case '&':
						case '^':
							buf.append(ch);
							accepted = true;
							lastOper = ch;
							state = state_oper_dbl;
							break;
						case '>':
						case '<':
						case '=':
							buf.append(ch);
							accepted = true;
							lastOper = ch;
							state = state_oper_compeql;
					}
					break;
				case state_oper_dbl:
					if (ch == lastOper) {
						buf.append(ch);
						accepted = true;
						if (ch == '>') {
							state = state_oper_tpl;
						} else {
							tokens.add(buf.toString());
							buf.delete(0, buf.length());
							states.add(state_lookup[state]);
							state = state_clear;
							lastOper = 0;
						}
					}
					break;
				case state_oper_tpl:
					if (ch == lastOper) {
						buf.append(ch);
						accepted = true;
						tokens.add(buf.toString());
						buf.delete(0, buf.length());
						states.add(state_lookup[state]);
						state = state_clear;
						lastOper = 0;
					}
					break;
				case state_oper_compeql:
					if (ch == '=') {
						buf.append(ch);
						accepted = true;
						tokens.add(buf.toString());
						buf.delete(0, buf.length());
						states.add(state_lookup[state]);
						state = state_clear;
						lastOper = 0;
					}
					break;
				case state_string:
					buf.append(ch);
					accepted = true;
					if (ch == '\'') {
						tokens.add(buf.toString());
						buf.delete(0, buf.length());
						states.add(state_lookup[state]);
						state = state_clear;
					}
					break;
				case state_string2:
					buf.append(ch);
					accepted = true;
					if (ch == '\"') {
						tokens.add(buf.toString());
						buf.delete(0, buf.length());
						states.add(state_lookup[state]);
						state = state_clear;
					}
					break;
				default: // case state_clear
					if (('0' <= ch && ch <= '9') || (ch == '.')) {
						state = state_num;
					} else if (('a' <= ch && ch <= 'z') || ('A' <= ch && ch <= 'Z')) {
						state = state_ident;
					} else if (ch == '\'') {
						buf.append(ch);
						state = state_string;
						continue;
					} else if (ch == '\"') {
						buf.append(ch);
						state = state_string2;
						continue;
					} else {
						switch (ch) {
							case '+':
							case '-':
							case '*':
							case '/':
							case '%':
							case '~':
							case '!':
							case ',':
							case '(':
							case ')':
							case '[':
							case ']':
							case '{':
							case '}':
							case '|':
							case '&':
							case '^':
							case '>':
							case '<':
							case '=':
							case '?':
							case ':':
								state = state_oper;
								lastOper = 0;
								break;
							default:
								throw new RuntimeException();
						}
					}
					i--;
					continue;
			}

			if (!accepted) {
				i--;
				if (buf.length() > 0) {
					tokens.add(buf.toString());
					buf.delete(0, buf.length());
					states.add(state_lookup[state]);
				}
				state = state_clear;
			}
		}
		if (buf.length() > 0) {
			tokens.add(buf.toString());
			buf.delete(0, buf.length());
			states.add(state_lookup[state]);
		}

		if (states.size() != tokens.size())
			throw new RuntimeException("Error in parsing expression " + line);

		int[] state_array = new int[states.size()];
		for (int i = 0; i < state_array.length; i++) {
			state_array[i] = states.get(i);
			String curTok = tokens.get(i);
			if (state_array[i] == state_oper && !acceptedOperatorSet.contains(curTok))
				throw new RuntimeException("Unacceptable operator " + curTok);
		}
		String[] token_array = tokens.toArray(new String[tokens.size()]);
		return new TokenManager(token_array, state_array);
	}

	public int getCurrentTokenType() {
		return hasMoreTokens() ? token_types[pos] : state_invalid;
	}

	public String getCurrentToken() {
		return hasMoreTokens() ? tokens[pos] : null;
	}

	public void consumeToken() {
		if (hasMoreTokens()) pos++;
	}

	public boolean hasMoreTokens() {
		return pos < token_types.length;
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder("[");
		if (tokens != null && tokens.length > 0) {
			buf.append(tokens[0]);
			for (int i = 1; i < tokens.length; i++)
				buf.append(", " + tokens[i]);
		}
		buf.append("]");
		return buf + "@" + pos;
	}
}

enum ASTNodeType {
	IDENTIFIER, LITERAL_INT, LITERAL_FLOAT, LITERAL_STRING, FUN_CALL,
	PAREN, CURLY_PAREN, UNARY, POWER, MULTIPLICATIVE, ADDITIVE, INSET,
	BITSHIFT, COMPARISON, BITWISE, BOOLEAN, QMARK_CHOICE;

	public boolean isLiteral() {
		return this == LITERAL_INT || this == LITERAL_FLOAT || this == LITERAL_STRING;
	}
}

enum DataType {
	UNKNOWN, NUMERIC, STRING, BOOLEAN, SET
}

class ASTNode {
	private ASTNodeType type;
	private DataType dataType = DataType.UNKNOWN;
	private String token;
	private List<ASTNode> children;

	public ASTNode(String tok, ASTNodeType t) {
		token = tok; type = t;
	}

	public ASTNode(String tok, ASTNodeType t, DataType dt) {
		token = tok; type = t; dataType = dt;
	}

	public void addChildren(ASTNode... nodes) {
		if (children == null)
			children = new ArrayList<ASTNode>();
		for (ASTNode node: nodes)
			children.add(node);
	}

	public void addChildren(Collection<ASTNode> nodes) {
		if (children == null)
			children = new ArrayList<ASTNode>();
		for (ASTNode node: nodes)
			children.add(node);
	}

	public String toString() {
		return children == null ? token : token + children;
	}

	public ASTNodeType getNodeType() {
		return type;
	}

	public DataType getDataType() {
		return dataType;
	}

	public void setDataType(DataType dt) {
		dataType = dt;
	}

	public void pushDataType(DataType dt) {
		if (dataType == DataType.UNKNOWN)
			dataType = dt;
		if (children != null)
			for (ASTNode child: children)
				child.pushDataType(dt);
	}

	public String getToken() {
		return token;
	}

	public ASTNode getChild(int i) {
		return children == null ? null : children.get(i);
	}

	public int getNumChildren() {
		return children == null ? 0 : children.size();
	}

	public List<ASTNode> getChildren() {
		return children;
	}
}

class ASTBuilder {
	public static final ASTNode parseExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseBooleanExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || !tok.equals("?"))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild1 = parseBooleanExpr(tokmgr);
			String tok2 = tokmgr.getCurrentToken();
			if (tok2 == null || !tok2.equals(":"))
				throw new RuntimeException("Missing ':' after '?'");
			tokmgr.consumeToken();
			ASTNode rightChild2 = parseBooleanExpr(tokmgr);
			ASTNode node = new ASTNode(tok2, ASTNodeType.QMARK_CHOICE);
			node.addChildren(leftChild, rightChild1, rightChild2);
			leftChild = node;
		} while (true);
		return leftChild;
	}

	static final ASTNode parseBooleanExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseBitwiseExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("&&") && !tok.equals("||") && !tok.equals("^^")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseBitwiseExpr(tokmgr);
			ASTNode node = new ASTNode(tok, ASTNodeType.BOOLEAN, DataType.BOOLEAN);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while (true);
		return leftChild;
	}

	static final ASTNode parseBitwiseExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseComparisonExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("&") && !tok.equals("|") && !tok.equals("^")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseComparisonExpr(tokmgr);
			ASTNode node = new ASTNode(tok, ASTNodeType.BITWISE);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while (true);
		return leftChild;
	}

	static final ASTNode parseComparisonExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseBitShiftExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("==") && !tok.equals("<") && !tok.equals("<=") && !tok.equals(">") && !tok.equals(">=")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseBitShiftExpr(tokmgr);
			ASTNode node = new ASTNode(tok, ASTNodeType.COMPARISON, DataType.BOOLEAN);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while (true);
		return leftChild;
	}

	static final ASTNode parseBitShiftExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseInSetExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("<<") && !tok.equals(">>") && !tok.equals(">>>")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseInSetExpr(tokmgr);
			ASTNodeType
				leftChildNodeType = leftChild.getNodeType(),
				rightChildNodeType = rightChild.getNodeType();
			if (leftChildNodeType.isLiteral() && rightChildNodeType.isLiteral()) {
				if (leftChildNodeType != ASTNodeType.LITERAL_INT && rightChildNodeType != ASTNodeType.LITERAL_INT)
					throw new RuntimeException("Bitshift expression must be executed on integer variables");
				long
					left = Long.parseLong(leftChild.getToken()),
					right = Long.parseLong(rightChild.getToken());
				if (tok.equals("<<"))
					right = left << right;
				else if (tok.equals(">>"))
					right = left >> right;
				else if (tok.equals(">>>"))
					right = left >>> right;
				return new ASTNode(String.valueOf(right), ASTNodeType.LITERAL_INT, DataType.NUMERIC);
			}
			ASTNode node = new ASTNode(tok, ASTNodeType.BITSHIFT, DataType.NUMERIC);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while (true);
		return leftChild;
	}

	static final ASTNode parseInSetExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseAdditiveExpr(tokmgr);
		String tok = tokmgr.getCurrentToken();
		if (tok != null && tok.equals("%%")) {
			tokmgr.consumeToken();
			ASTNode rightChild = parseCurlyParenExpr(tokmgr);
			ASTNode node = new ASTNode(tok, ASTNodeType.INSET, DataType.BOOLEAN);
			node.addChildren(leftChild, rightChild);
			return node;
		}
		return leftChild;
	}

	static final ASTNode parseAdditiveExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseMultiplicativeExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("+") && !tok.equals("-")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseMultiplicativeExpr(tokmgr);
			ASTNodeType
				leftChildNodeType = leftChild.getNodeType(),
				rightChildNodeType = rightChild.getNodeType();
			if (leftChildNodeType.isLiteral() && rightChildNodeType.isLiteral()) {
				if (leftChildNodeType == ASTNodeType.LITERAL_STRING || rightChildNodeType == ASTNodeType.LITERAL_STRING)
					throw new RuntimeException("Addition or subtraction expression must be executed on numerical variables");
				double
					left = Double.parseDouble(leftChild.getToken()),
					right = Double.parseDouble(rightChild.getToken());
				if (tok.equals("+"))
					right = left + right;
				else if (tok.equals("-"))
					right = left - right;
				return new ASTNode(String.valueOf(right), ASTNodeType.LITERAL_FLOAT, DataType.NUMERIC);
			}
			ASTNode node = new ASTNode(tok, ASTNodeType.ADDITIVE);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while(true);
		return leftChild;
	}

	static final ASTNode parseMultiplicativeExpr(TokenManager tokmgr) {
		ASTNode leftChild = parsePowerExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || (!tok.equals("*") && !tok.equals("/") && !tok.equals("%")))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parsePowerExpr(tokmgr);
			ASTNodeType
				leftChildNodeType = leftChild.getNodeType(),
				rightChildNodeType = rightChild.getNodeType();
			if (leftChildNodeType.isLiteral() && rightChildNodeType.isLiteral()) {
				if (leftChildNodeType == ASTNodeType.LITERAL_STRING || rightChildNodeType == ASTNodeType.LITERAL_STRING)
					throw new RuntimeException("Multiplicative, division, or remainder expression must be executed on numerical variables");
				double
					left = Double.parseDouble(leftChild.getToken()),
					right = Double.parseDouble(rightChild.getToken());
				if (tok.equals("*"))
					right = left * right;
				else if (tok.equals("/"))
					right = left / right;
				else if (tok.equals("%"))
					right = left % right;
				return new ASTNode(String.valueOf(right), ASTNodeType.LITERAL_FLOAT, DataType.NUMERIC);
			}
			ASTNode node = new ASTNode(tok, ASTNodeType.MULTIPLICATIVE, DataType.NUMERIC);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while(true);
		return leftChild;
	}

	static final ASTNode parsePowerExpr(TokenManager tokmgr) {
		ASTNode leftChild = parseUnaryExpr(tokmgr);
		do {
			String tok = tokmgr.getCurrentToken();
			if (tok == null || !tok.equals("**"))
				break;
			tokmgr.consumeToken();
			ASTNode rightChild = parseUnaryExpr(tokmgr);
			ASTNodeType
				leftChildNodeType = leftChild.getNodeType(),
				rightChildNodeType = rightChild.getNodeType();
			if (leftChildNodeType.isLiteral() && rightChildNodeType.isLiteral()) {
				if (leftChildNodeType == ASTNodeType.LITERAL_STRING || rightChildNodeType == ASTNodeType.LITERAL_STRING)
					throw new RuntimeException("Power expression must be executed on numerical variables");
				double
					left = Double.parseDouble(leftChild.getToken()),
					right = Double.parseDouble(rightChild.getToken());
				return new ASTNode(String.valueOf(pow(left, right)), ASTNodeType.LITERAL_FLOAT, DataType.NUMERIC);
			}
			ASTNode node = new ASTNode(tok, ASTNodeType.POWER, DataType.NUMERIC);
			node.addChildren(leftChild, rightChild);
			leftChild = node;
		} while(true);
		return leftChild;
	}

	static final ASTNode parseUnaryExpr(TokenManager tokmgr) {
		String tok = tokmgr.getCurrentToken();
		if (tok != null && (tok.equals("+") || tok.equals("-") || tok.equals("!") || tok.equals("~"))) {
			tokmgr.consumeToken();
			ASTNode node = new ASTNode(tok, ASTNodeType.UNARY);
			ASTNode child = parseCurlyParenExpr(tokmgr);
			if (child == null)
				throw new RuntimeException("Expecting an expression!");
			ASTNodeType childNodeType = child.getNodeType();
			if (childNodeType.isLiteral()) {
				if (childNodeType == ASTNodeType.LITERAL_INT || childNodeType == ASTNodeType.LITERAL_FLOAT) {
					node = new ASTNode(tok + child.getToken(), childNodeType, child.getDataType());
				} else
					throw new RuntimeException("Not a valid expression: " + tok + child.getToken());
			} else {
				node.setDataType(child.getDataType());
				node.addChildren(child);
			}
			return node;
		}
		return parseCurlyParenExpr(tokmgr);
	}

	static final ASTNode parseCurlyParenExpr(TokenManager tokmgr) {
		String tok = tokmgr.getCurrentToken();
		if (tok != null && tok.equals("{")) {
			tokmgr.consumeToken();
			ASTNode node = new ASTNode(tok, ASTNodeType.CURLY_PAREN, DataType.SET);
			DataType dt = null;
			while (true) {
				ASTNode curNode = parseExpr(tokmgr);
				if (curNode == null)
					break;
				if (!curNode.getNodeType().isLiteral())
					throw new RuntimeException("Expressions in curly parentheses must be literals!");
				node.addChildren(curNode);
				if (dt == null)
					dt = curNode.getDataType();
				else
					if (dt != curNode.getDataType())
						throw new RuntimeException("Literals in curly parentheses must be of the same type!");
				tok = tokmgr.getCurrentToken();
				if (!tok.equals(","))
					break;
				tokmgr.consumeToken();
			}
			tok = tokmgr.getCurrentToken();
			if (!tok.equals("}"))
				throw new RuntimeException("Expecting right curly parenthesis!");
			tokmgr.consumeToken();
			return node;
		}
		return parseParenExpr(tokmgr);
	}

	static final ASTNode parseParenExpr(TokenManager tokmgr) {
		String tok = tokmgr.getCurrentToken();
		if (tok != null && tok.equals("(")) {
			tokmgr.consumeToken();
			ASTNode node = new ASTNode(tok, ASTNodeType.PAREN);
			while (true) {
				ASTNode curNode = parseExpr(tokmgr);
				if (curNode == null)
					break;
				node.addChildren(curNode);
				tok = tokmgr.getCurrentToken();
				if (!tok.equals(","))
					break;
				tokmgr.consumeToken();
			}
			tok = tokmgr.getCurrentToken();
			if (!tok.equals(")"))
				throw new RuntimeException("Expecting right parenthesis!");
			tokmgr.consumeToken();
			return node;
		}
		return parseSingularExpr(tokmgr);
	}

	static final ASTNode parseSingularExpr(TokenManager tokmgr) {
		String tok = tokmgr.getCurrentToken();
		if (tok == null)
			return null;
		int curTokenType = tokmgr.getCurrentTokenType();
		switch(curTokenType) {
			case TokenManager.state_long:
				tokmgr.consumeToken();
				return new ASTNode(tok, ASTNodeType.LITERAL_INT, DataType.NUMERIC);
			case TokenManager.state_float:
				tokmgr.consumeToken();
				return new ASTNode(tok, ASTNodeType.LITERAL_FLOAT, DataType.NUMERIC);
			case TokenManager.state_ident:
				tokmgr.consumeToken();
				String tok2 = tokmgr.getCurrentToken();
				if (tok2 != null && tok2.equals("(")) {
					ASTNode paren = parseParenExpr(tokmgr);
					ASTNode node = new ASTNode(tok, ASTNodeType.FUN_CALL);
					node.addChildren(paren.getChildren());
					return node;
				}
				return new ASTNode(tok, ASTNodeType.IDENTIFIER);
			case TokenManager.state_string:
				tokmgr.consumeToken();
				tok = "\"" + tok.substring(1, tok.length()-1) + "\"";
				return new ASTNode(tok, ASTNodeType.LITERAL_STRING, DataType.STRING);
			case TokenManager.state_invalid:
				throw new RuntimeException("Unexpected end of expression!");
			default:
				return null;
		}
	}

	private static final void getTerminalNodes(ASTNode node, List<ASTNode> termNodes) {
		List<ASTNode> children = node.getChildren();
		if (children == null)
			termNodes.add(node);
		else {
			for (ASTNode child: children)
				getTerminalNodes(child, termNodes);
		}
	}

	public static final List<ASTNode> getTerminalNodes(ASTNode node) {
		List<ASTNode> termNodes = new ArrayList<ASTNode>();
		getTerminalNodes(node, termNodes);
		return termNodes;
	}

	public static final void inferTypes(ASTNode node) {
		List<ASTNode> children = node.getChildren();
		if (children != null)
			for (ASTNode child: children)
				inferTypes(child);
		ASTNode left, right;
		DataType leftDt, rightDt;
		switch (node.getNodeType()) {
			case ADDITIVE:
				left = children.get(0);
				right = children.get(1);
				leftDt = left.getDataType();
				rightDt = right.getDataType();
				if (leftDt == DataType.UNKNOWN) {
					if (rightDt == DataType.UNKNOWN) {
						left.pushDataType(DataType.NUMERIC);
						right.pushDataType(DataType.NUMERIC);
						node.setDataType(DataType.NUMERIC);
						System.out.println("Cannot infer types of the expression " + node + ", assuming numeric.");
					} else {
						if (rightDt != DataType.NUMERIC || rightDt != DataType.STRING)
							throw new RuntimeException("Incompatible operation for " + node);
						left.pushDataType(rightDt);
						node.setDataType(rightDt);
					}
				} else {
					if (rightDt == DataType.UNKNOWN) {
						if (leftDt != DataType.NUMERIC || leftDt != DataType.STRING)
							throw new RuntimeException("Incompatible operation for " + node);
						right.pushDataType(leftDt);
						node.setDataType(leftDt);
					} else {
						if (leftDt != rightDt)
							throw new RuntimeException("Conflicting data types of the expression " + node);
						String tok = node.getToken();
						if (leftDt == DataType.STRING && !tok.equals("+"))
							throw new RuntimeException("Subtraction cannot be used in strings " + node);
						else if (leftDt != DataType.NUMERIC)
							throw new RuntimeException("Incompatible operation for " + node);
					}
				}
				break;
			case BITSHIFT:
			case BITWISE:
			case MULTIPLICATIVE:
			case POWER:
				left = children.get(0);
				right = children.get(1);
				leftDt = left.getDataType();
				rightDt = right.getDataType();
				if (leftDt == DataType.UNKNOWN) {
					if (rightDt == DataType.UNKNOWN) {
						throw new RuntimeException("Cannot infer types of the expression " + node);
					} else {
						if (rightDt != DataType.NUMERIC)
							throw new RuntimeException("Incompatible operation for " + node);
						left.pushDataType(rightDt);
						node.setDataType(rightDt);
					}
				} else {
					if (rightDt == DataType.UNKNOWN) {
						if (leftDt != DataType.NUMERIC)
							throw new RuntimeException("Incompatible operation for " + node);
						right.pushDataType(leftDt);
						node.setDataType(leftDt);
					} else {
						if (leftDt != rightDt)
							throw new RuntimeException("Conflicting data types of the expression " + node);
						if (leftDt != DataType.NUMERIC)
							throw new RuntimeException("Incompatible operation for " + node);
					}
				}
				break;
			case COMPARISON:
				left = children.get(0);
				right = children.get(1);
				leftDt = left.getDataType();
				rightDt = right.getDataType();
				if (leftDt == DataType.UNKNOWN) {
					if (rightDt == DataType.UNKNOWN) {
						if (node.getToken().equals("==")) {
							left.pushDataType(DataType.STRING);
							right.pushDataType(DataType.STRING);
							node.setDataType(DataType.BOOLEAN);
							System.out.println("Cannot infer types of the expression " + node + ", assuming string.");
						} else {
							left.pushDataType(DataType.NUMERIC);
							right.pushDataType(DataType.NUMERIC);
							node.setDataType(DataType.BOOLEAN);
							System.out.println("Cannot infer types of the expression " + node + ", assuming numeric.");
						}
					} else {
						left.pushDataType(rightDt);
						node.setDataType(DataType.BOOLEAN);
					}
				} else {
					if (rightDt == DataType.UNKNOWN) {
						right.pushDataType(leftDt);
						node.setDataType(DataType.BOOLEAN);
					} else {
						if (leftDt != rightDt)
							throw new RuntimeException("Conflicting data types of the expression " + node);
						if (leftDt != DataType.NUMERIC)
							throw new RuntimeException("Incompatible operation for " + node);
					}
				}
				break;
			case BOOLEAN:
				left = children.get(0);
				right = children.get(1);
				leftDt = left.getDataType();
				rightDt = right.getDataType();
				if (leftDt == DataType.UNKNOWN) {
					if (rightDt == DataType.UNKNOWN) {
						throw new RuntimeException("Cannot infer types of the expression " + node);
					} else {
						if (rightDt != DataType.BOOLEAN)
							throw new RuntimeException("Incompatible operation for " + node);
						left.pushDataType(rightDt);
						node.setDataType(rightDt);
					}
				} else {
					if (rightDt == DataType.UNKNOWN) {
						if (leftDt != DataType.BOOLEAN)
							throw new RuntimeException("Incompatible operation for " + node);
						right.pushDataType(leftDt);
						node.setDataType(leftDt);
					} else {
						if (leftDt != rightDt)
							throw new RuntimeException("Conflicting data types of the expression " + node);
						if (leftDt != DataType.BOOLEAN)
							throw new RuntimeException("Incompatible operation for " + node);
					}
				}
				break;
			case FUN_CALL:
				String tok = node.getToken();
				if (TokenManager.numericFunCall.contains(tok))
					node.pushDataType(DataType.NUMERIC);
				else if (TokenManager.stringFunCall.contains(tok))
					node.pushDataType(DataType.STRING);
				else if (TokenManager.stringReturningBoolFunCall.contains(tok)) {
					node.pushDataType(DataType.STRING);
					node.setDataType(DataType.BOOLEAN);
				} else if (TokenManager.stringReturningIntFunCall.contains(tok)) {
					node.pushDataType(DataType.STRING);
					node.setDataType(DataType.NUMERIC);
				} else if (tok.equals("split")) {
					node.setDataType(DataType.SET);
					node.getChild(0).setDataType(DataType.STRING);
				} else if (tok.equals("file"))
					node.setDataType(DataType.SET);
				break;
			case UNARY:
			case PAREN:
				left = children.get(0);
				leftDt = left.getDataType();
				node.setDataType(leftDt);
				break;
			case QMARK_CHOICE:
				left = children.get(0);
				right = children.get(1);
				ASTNode right2 = children.get(2);
				leftDt = left.getDataType();
				rightDt = right.getDataType();
				DataType rightDt2 = right2.getDataType();

				if (rightDt == DataType.UNKNOWN) {
					if (rightDt2 == DataType.UNKNOWN) {
						throw new RuntimeException("Cannot infer types of the expression " + node);
					} else {
						right.pushDataType(rightDt2);
						node.setDataType(rightDt2);
					}
				} else {
					if (rightDt2 == DataType.UNKNOWN) {
						right2.pushDataType(rightDt);
						node.setDataType(rightDt);
					} else {
						if (rightDt != rightDt2)
							throw new RuntimeException("Conflicting data types of the expression " + node);
						node.setDataType(rightDt);
					}
				}
				break;
			default:
		}
	}
}

class JavaProgram {
	Set<String>
		importSet = new HashSet<String>(),
		importStaticSet = new HashSet<String>();
	List<String> fields = new ArrayList<String>();
	StringBuilder
		constructor = new StringBuilder(),
		method = new StringBuilder();
	boolean needLoadMethod = false;
}

class CodeGenerator {
	static int count = 0;
	static final String
		className = "MyFilter",
		methodName = "accept",
		objPrefix = "_obj",
		ln = System.getProperty("line.separator"),
		tab = "    ";

	public static final String generate(ASTNode node, IdentityHashMap<ASTNode, String> symbolTable) {
		JavaProgram prog = new JavaProgram();
		prog.importStaticSet.add("java.lang.Math");
		prog.importSet.add("java.util.*");
		generateConstructor(node, prog, symbolTable);
		generateMethod(node, prog, symbolTable);

		if (prog.needLoadMethod)
			prog.importSet.add("com.csvreader.CsvReader");

		StringBuilder buf = new StringBuilder();
		for (String imp: prog.importSet) {
			buf.append("import " + imp + ";" + ln);
		}
		for (String imp: prog.importStaticSet) {
			buf.append("import static " + imp + ".*;" + ln);
		}
		buf.append("public class "+ className + " {" + ln);
		for (String field: prog.fields) {
			buf.append(tab + "static " + field + ln);
		}
		// Add constructor code, if any
		if (prog.constructor.length() > 0) {
			buf.append(tab + "static {" + ln);
			buf.append(prog.constructor);
			buf.append(tab + "}" + ln);
		}
		buf.append(tab + "public static final boolean " + methodName + "(String[] tokens) {" + ln);
		// Add method code
		buf.append(tab + tab + "return " + prog.method + ";" + ln);
		buf.append(tab + "}" + ln);

		buf.append(tab + "public static final <T> Set<T> makeSet(T[] tokens) {" + ln);
		buf.append(tab + tab + "Set<T> set = new HashSet<T>();" + ln);
		buf.append(tab + tab + "for(T t: tokens) set.add(t);" + ln);
		buf.append(tab + tab + "return set;" + ln);
		buf.append(tab + "}" + ln);

		buf.append(tab + "public static final <T> boolean nonEmptyIntersect(Set<T> s1, Set<T> s2) {" + ln);
		buf.append(tab + tab + "Set<T> set = new HashSet<T>(); set.addAll(s1);" + ln);
		buf.append(tab + tab + "set.retainAll(s2);" + ln);
		buf.append(tab + tab + "return set.size() > 0;" + ln);
		buf.append(tab + "}" + ln);

		if (prog.needLoadMethod) {
			buf.append(tab + "public static final Set<String> load(String filename, int colNo, char sep) {" + ln);
			buf.append(tab + tab + "Set<String> set = new HashSet<String>();" + ln);
			buf.append(tab + tab + "CsvReader reader = null;" + ln);
			buf.append(tab + tab + "try {" + ln);
			buf.append(tab + tab + tab + "reader = new CsvReader(filename, new char[] {sep});" + ln);
			buf.append(tab + tab + tab + "while(reader.readRecord()) {" + ln);
			buf.append(tab + tab + tab + tab + "set.add(reader.getValues()[colNo]);" + ln);
			buf.append(tab + tab + tab + "}" + ln);
			buf.append(tab + tab + "} catch (Exception e) {" + ln);
			buf.append(tab + tab + tab + "e.printStackTrace();" + ln);
			buf.append(tab + tab + "} finally {" + ln);
			buf.append(tab + tab + tab + "if (reader != null)	try { reader.close(); } catch (Exception e) {}" + ln);
			buf.append(tab + tab + "}" + ln);
			buf.append(tab + tab + "return set;" + ln);
			buf.append(tab + "}" + ln);
		}
		buf.append("}" + ln);
		return buf.toString();
	}

	static final void generateConstructor(ASTNode node, JavaProgram prog, IdentityHashMap<ASTNode, String> symbolTable) {
		String objName, tok;
		List<ASTNode> children = node.getChildren();
		if (children != null)
			for (ASTNode child: children)
				generateConstructor(child, prog, symbolTable);
		switch (node.getNodeType()) {
			case CURLY_PAREN:
				objName = objPrefix+count;
				symbolTable.put(node, objName);
				count++;
				prog.importSet.add("java.util.Set");
				prog.importSet.add("java.util.HashSet");
				switch (children.get(0).getDataType()) {
					case NUMERIC:
						prog.fields.add("Set<Double> " + objName + " = new HashSet<Double>();");
						break;
					default:
						prog.fields.add("Set<String> " + objName + " = new HashSet<String>();");
						break;
				}
				for (ASTNode child: children) {
					if (child.getNodeType().isLiteral()) {
						prog.constructor.append(tab + tab + objName + ".add(" + child.getToken() + ");" + ln);
					} else {
						prog.constructor.append(tab + tab + objName + ".addAll(" + symbolTable.get(child) + ");" + ln);
					}
				}
				break;
			case FUN_CALL:
				tok = node.getToken();
				if (!tok.equals("file"))
					break;
				String fn = children.get(0).getToken();
				int colNo = children.size() > 1 ? Integer.parseInt(children.get(1).getToken()) : 0;
				char sep = children.size() > 2 ? children.get(2).getToken().charAt(0) : ',';
				objName = objPrefix+count;
				symbolTable.put(node, objName);
				count++;
				prog.importSet.add("java.util.Set");
				prog.importSet.add("java.util.HashSet");
				prog.needLoadMethod = true;
				prog.fields.add("Set<String> " + objName + ";");
				prog.constructor.append(tab + tab + objName + " = load(" + fn + ", " + colNo + ", '" + sep + "');" + ln);
				break;
			default:
		}
	}

	static final void generateMethod(ASTNode node, JavaProgram prog, IdentityHashMap<ASTNode, String> symbolTable) {
		generateMethod(node, prog, symbolTable, "tokens");
	}

	static final void generateMethod(ASTNode node, JavaProgram prog, IdentityHashMap<ASTNode, String> symbolTable,
		String tokenKeyword) {
		String objName, tok;
		List<ASTNode> children = node.getChildren();
		switch (node.getNodeType()) {
			case POWER:
				prog.method.append("pow(");
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				prog.method.append(",");
				generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
				prog.method.append(")");
				break;
			case INSET: // TODO
				objName = symbolTable.get(children.get(1));
				if (objName == null) {
					if (children.get(1).getDataType() != DataType.SET)
						throw new RuntimeException();
					prog.method.append("nonEmptyIntersect(");
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append(",");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
					prog.method.append(")");
				} else {
					prog.method.append(objName + ".contains(");
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append(")");
				}
				break;
			case UNARY:
				prog.method.append(node.getToken());
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				break;
			case PAREN:
				prog.method.append("(");
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				prog.method.append(")");
				break;
			case IDENTIFIER:
				objName = symbolTable.get(node);
				if (node.getDataType() == DataType.NUMERIC)
					prog.method.append("Double.parseDouble(");
				prog.method.append(tokenKeyword + "[" + objName + "]");
				if (node.getDataType() == DataType.NUMERIC)
					prog.method.append(")");
				break;
			case CURLY_PAREN:
				objName = symbolTable.get(node);
				prog.method.append(objName);
				break;
			case FUN_CALL:
				tok = node.getToken();
				if (tok.equals("file")) {
					objName = symbolTable.get(node);
					prog.method.append(objName);
					break;
				}
				if (tok.equals("substring") || tok.equals("startsWith") || tok.equals("replaceAll") || tok.equals("replaceFirst")
						|| tok.equals("indexOf") || tok.equals("lastIndexOf")) {
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append("." + tok + "(");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
					if (children.size() > 2) {
						prog.method.append(",");
						generateMethod(children.get(2), prog, symbolTable, tokenKeyword);
					}
					prog.method.append(")");
					break;
				}
				if (tok.equals("endsWith") || tok.equals("matches") || tok.equals("indexOf") || tok.equals("lastIndexOf") ||
					tok.equals("contains")) {
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append("." + tok + "(");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
					prog.method.append(")");
					break;
				}
				if (tok.equals("upper")) tok = "toUpperCase";
				if (tok.equals("lower")) tok = "toLowerCase";
				if (tok.equals("toUpperCase") || tok.equals("toLowerCase") || tok.equals("trim")) {
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append("." + tok + "()");
					break;
				}
				if (tok.equals("split")) { //TODO
					prog.method.append("makeSet(");
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append("." + tok + "(");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
					if (children.size() > 2) {
						prog.method.append(",");
						generateMethod(children.get(2), prog, symbolTable, tokenKeyword);
					}
					prog.method.append("))");
					break;
				}
				prog.method.append(tok);
				prog.method.append("(");
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				for (int i = 1; i < children.size(); i++) {
					prog.method.append(",");
					generateMethod(children.get(i), prog, symbolTable, tokenKeyword);
				}
				prog.method.append(")");
				break;
			case COMPARISON:
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				boolean isStr = children.get(0).getDataType() == DataType.STRING;
				if (isStr) {
					prog.method.append(".equals(");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
					prog.method.append(")");
				} else {
					prog.method.append(" " + node.getToken() + " ");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
				}
				break;
			case QMARK_CHOICE:
				generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
				prog.method.append(" ? ");
				generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
				prog.method.append(" : ");
				generateMethod(children.get(2), prog, symbolTable, tokenKeyword);
				break;
			default:
				if (children != null) {
					generateMethod(children.get(0), prog, symbolTable, tokenKeyword);
					prog.method.append(" " + node.getToken() + " ");
					generateMethod(children.get(1), prog, symbolTable, tokenKeyword);
				} else
					prog.method.append(node.getToken());
		}
	}
}

// Code for on-the-fly compilation

class Source extends SimpleJavaFileObject {
	private final String content;

	Source(String name, Kind kind, String content) {
		super(URI.create("memo:///" + name.replace('.', '/') + kind.extension), kind);
		this.content = content;
	}

	@Override
	public CharSequence getCharContent(boolean ignore) {
		return this.content;
	}
}

class Output extends SimpleJavaFileObject {
	private final ByteArrayOutputStream baos = new ByteArrayOutputStream();

	Output(String name, Kind kind) {
		super(URI.create("memo:///" + name.replace('.', '/') + kind.extension), kind);
	}

	byte[] toByteArray() {
		return this.baos.toByteArray();
	}

	@Override
	public ByteArrayOutputStream openOutputStream() {
		return this.baos;
	}
}

class MemoryFileManager extends ForwardingJavaFileManager<JavaFileManager> {
	final Map<String, Output> map = new HashMap<String, Output>();

	MemoryFileManager(JavaCompiler compiler) {
		super(compiler.getStandardFileManager(null, null, null));
	}

	@Override
	public Output getJavaFileForOutput
	(Location location, String name, Kind kind, FileObject source) {
		Output mc = new Output(name, kind);
		this.map.put(name, mc);
		return mc;
	}
}

class MemoryClassLoader extends ClassLoader {
	private final JavaCompiler compiler = ToolProvider.getSystemJavaCompiler();
	private final MemoryFileManager manager = new MemoryFileManager(this.compiler);

	public MemoryClassLoader(String classname, String filecontent) {
		this(Collections.singletonMap(classname, filecontent));
	}

	public MemoryClassLoader(Map<String, String> map) {
		List<Source> list = new ArrayList<Source>();
		for (Map.Entry<String, String> entry : map.entrySet()) {
			list.add(new Source(entry.getKey(), Kind.SOURCE, entry.getValue()));
		}
		List<String> optionList = new ArrayList<String>();
		optionList.addAll(Arrays.asList("-cp",System.getProperty("java.class.path")));
		this.compiler.getTask(null, this.manager, null, optionList, null, list).call();
	}

	@Override
	protected Class<?> findClass(String name) throws ClassNotFoundException {
		synchronized (this.manager) {
			Output mc = this.manager.map.remove(name);
			if (mc != null) {
				byte[] array = mc.toByteArray();
				return defineClass(name, array, 0, array.length);
			}
		}
		return super.findClass(name);
	}
}
