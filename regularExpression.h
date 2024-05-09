#ifndef REGULAREXPRESSION_H
#define REGULAREXPRESSION_H

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <stack>
#include <utility>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include "constants.h"

#define NDEBUG

struct SymbolOrEpsilon
{
	Symbol c;
	constexpr SymbolOrEpsilon(): c(Constants::Epsilon) {}
	constexpr SymbolOrEpsilon(Symbol c) : c(c) {}
	constexpr SymbolOrEpsilon(std::string_view s) : c(s[0]) {}
	static constexpr SymbolOrEpsilon epsilon() noexcept { return {}; }
	static constexpr bool isBegin(Symbol c) noexcept
	{
		return c == Constants::Epsilon || !Constants::isSpecial(c);
	}
	static constexpr bool isEnd(Symbol c) noexcept
	{
		return isBegin(c);
	}
	auto operator<=>(Symbol c) const noexcept
	{
		return this->c <=> c;
	}
	bool operator==(Symbol c) const noexcept
	{
		return this->c == c;
	}
	operator Symbol() const noexcept
	{
		return c;
	}
	auto operator<=>(const SymbolOrEpsilon&) const noexcept = default;
	bool operator==(const SymbolOrEpsilon&) const noexcept = default;
};

struct SymbolPair
{
	SymbolOrEpsilon coords[2];
	static constexpr SymbolPair epsilon() { return {}; }
	auto operator<=>(const SymbolPair&) const noexcept = default;
	bool operator==(const SymbolPair&) const noexcept = default;
};

struct WordPair
{
	Word coords[2];
	constexpr WordPair(): coords{Word(1, Constants::Epsilon), Word(1, Constants::Epsilon)} {}
	constexpr WordPair(std::string_view wp)
	{
		wp.remove_prefix(1); // Constants::BaseElementBegin
		wp.remove_suffix(1); // Constants::BaseElementEnd
		std::size_t delimPos = wp.find(Constants::BaseElementDelim);
		if(delimPos == std::string_view::npos)
			throw std::runtime_error("Bad regular expression: missing delimiter in base element");
		coords[0] = wp.substr(0, delimPos);
		coords[1] = wp.substr(delimPos + 1);
		if(coords[0].empty() || coords[1].empty() || coords[1].find(Constants::BaseElementDelim) != std::string::npos)
			throw std::runtime_error("Bad regular expression: missing coordinate in base element or there are extra delimiters");
	}
	static constexpr WordPair epsilon() { return {}; }
	static constexpr bool isBegin(Symbol c) noexcept
	{
		return c == Constants::BaseElementBegin;
	}
	static constexpr bool isEnd(Symbol c) noexcept
	{
		return c == Constants::BaseElementEnd;
	}
	auto operator<=>(const WordPair&) const noexcept = default;
	bool operator==(const WordPair&) const noexcept = default;

	SymbolPair operator()(std::size_t ind) const noexcept
	{
		auto proj = [](const Word& w, std::size_t ind) -> Symbol { return ind < w.size() ? w[ind] : Constants::Epsilon; };
		return {proj(coords[0], ind), proj(coords[1], ind)};
	}
};

struct Symbol_Word
{
	Symbol first;
	Word second;
	constexpr Symbol_Word(): first(Constants::Epsilon), second(1, Constants::Epsilon) {}
	Symbol_Word(Symbol first, const Word& second): first(first), second(second) {}
	Symbol_Word(Symbol first, Word&& second): first(first), second(std::move(second)) {}
	static constexpr Symbol_Word epsilon() { return {}; }
	auto operator<=>(const Symbol_Word&) const noexcept = default;
	bool operator==(const Symbol_Word&) const noexcept = default;
};

inline std::ostream& operator<<(std::ostream& os, const SymbolPair& sp) // [a,b]
{
	return os << Constants::BaseElementBegin << sp.coords[0] << Constants::BaseElementDelim << sp.coords[1] << Constants::BaseElementEnd;
}
inline std::ostream& operator<<(std::ostream& os, const Symbol_Word& sw) // [a,"b"]
{
	return os << Constants::BaseElementBegin << sw.first << Constants::BaseElementDelim << std::quoted(sw.second) << Constants::BaseElementEnd;
}
inline std::ostream & operator<<(std::ostream & os, const WordPair & wp) // ["a","b"]
{
	return os << Constants::BaseElementBegin << std::quoted(wp.coords[0]) << Constants::BaseElementDelim << std::quoted(wp.coords[1]) << Constants::BaseElementEnd;
}

inline std::istream& operator>>(std::istream& is, SymbolOrEpsilon& s) // a
{
	return is >> s.c;
}
inline std::istream& operator>>(std::istream& is, SymbolPair& sp) // [a,b]
{
	char ignore;
	return is >> ignore >> sp.coords[0] >> ignore >> sp.coords[1] >> ignore;
}
inline std::istream& operator>>(std::istream& is, Symbol_Word& sw)
{
	char ignore;
	return is >> ignore >> sw.first >> ignore >> std::quoted(sw.second) >> ignore;
}
inline std::istream& operator>>(std::istream& is, WordPair& wp) // ["a","b"]
{
	char ignore;
	return is >> ignore >> std::quoted(wp.coords[0]) >> ignore >> std::quoted(wp.coords[1]) >> ignore;
}

template<class BaseElement>
class RegularExpression
{
#ifndef NDEBUG
	std::string tokenizedRegex;
#endif
	std::string tokenizedRPN;
	std::vector<BaseElement> baseTokens;
	std::string tokenize(const std::string& regex)
	{
		std::string tokenizedRegex;
		tokenizedRegex.reserve(regex.size());
		for(std::size_t i = 0; i < regex.size(); i++)
		{
			if((BaseElement::isBegin(regex[i]) || regex[i] == Constants::EmptySet || regex[i] == Constants::OpenParenthesis) &&
			   i && regex[i - 1] != Constants::OpenParenthesis && regex[i - 1] != Constants::Union)
				tokenizedRegex.push_back(Constants::Concatenation);
			if(Constants::isOperator(regex[i]) || Constants::isParenthesis(regex[i]) || regex[i] == Constants::EmptySet)
				tokenizedRegex.push_back(regex[i]);
			else if(BaseElement::isBegin(regex[i]))
				if(auto baseStart = regex.begin() + i, baseEnd = std::find_if(baseStart, regex.end(), BaseElement::isEnd); baseEnd != regex.end())
				{
					baseTokens.emplace_back(std::string_view{baseStart, baseEnd + 1});
					tokenizedRegex.push_back(Constants::BasePlaceholder);
					i = std::distance(regex.begin(), baseEnd);
				}
				else throw std::runtime_error("Bad regular expression: unclosed base element");
			else throw std::runtime_error("Bad regular expression: unexpected character at position " + std::to_string(i));
		}
		return tokenizedRegex;
	}
	void produceRPN(const std::string& tokenizedRegex)
	{
		tokenizedRPN.reserve(tokenizedRegex.size());
		std::stack<Symbol> op;
		for(Symbol c : tokenizedRegex)
		{
			if(c == Constants::BasePlaceholder || c == Constants::EmptySet)
				tokenizedRPN.push_back(c);
			else if(Constants::isOperator(c))
			{
				while(!op.empty() && op.top() != Constants::OpenParenthesis && Constants::precedence(op.top()) >= Constants::precedence(c))
				{
					tokenizedRPN.push_back(op.top());
					op.pop();
				}
				op.push(c);
			}
			else if(c == Constants::OpenParenthesis) op.push(c);
			else if(c == Constants::CloseParenthesis)
			{
				while(!op.empty() && op.top() != Constants::OpenParenthesis)
				{
					tokenizedRPN.push_back(op.top());
					op.pop();
				}
				if(op.empty()) throw std::runtime_error("Bad regular expression: mismatched parentheses");
				op.pop();
			}
		}
		while(!op.empty())
		{
			if(op.top() == Constants::OpenParenthesis) throw std::runtime_error("Bad regular expression: mismatched parentheses");
			tokenizedRPN.push_back(op.top());
			op.pop();
		}
	}
public:
	RegularExpression(const std::string& regex)
	{
		if(regex.empty())
			throw std::runtime_error("Empty regular expression");
#	ifndef NDEBUG
		tokenizedRegex = tokenize(regex);
		produceRPN(tokenizedRegex);
#	else
		produceRPN(tokenize(regex));
#	endif
	}
#ifndef NDEBUG
	const std::string& TokenizedRegex() const noexcept { return tokenizedRegex; }
#endif
	const std::string& TokenizedReversePolishNotation() const noexcept { return tokenizedRPN; }
	const std::vector<BaseElement>& BaseTokens() const noexcept { return baseTokens; }
};

#endif
