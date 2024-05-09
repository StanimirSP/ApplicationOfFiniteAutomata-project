#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstddef>
#include <string>
#include <utility>
#include <cstdint>

//using State = std::size_t;
using State = std::uint32_t;
using Symbol = char;
using USymbol = unsigned char;
using Word = std::string;

namespace Constants
{
	struct AlwaysTrue
	{
		template<class... Args>
		constexpr bool operator()(Args&&...) const noexcept
		{
			return true;
		}
	};

	constexpr State InvalidState = -1;

	constexpr Symbol Epsilon = '_';
	constexpr Symbol BaseElementBegin = '[';
	constexpr Symbol BaseElementEnd = ']';
	constexpr Symbol BaseElementDelim = ',';
	constexpr Symbol Union = '|';
	constexpr Symbol Concatenation = '&';
	constexpr Symbol KleeneStar = '*';
	constexpr Symbol OpenParenthesis = '(';
	constexpr Symbol CloseParenthesis = ')';
	constexpr Symbol BasePlaceholder = '#';
	constexpr Symbol EmptySet = '@';
	constexpr Symbol ReplacementPos = '^';
	constexpr Symbol ReplacementStart = '<';
	constexpr Symbol ReplacementEnd = '>';
	static_assert(BaseElementBegin != BaseElementEnd, "BaseElementBegin and BaseElementEnd must be different");

	bool isSpecial(USymbol c);
	bool isOperator(USymbol c);
	bool isParenthesis(USymbol c);
	bool isForbidden(USymbol c);
	int precedence(Symbol c);
}

#endif
