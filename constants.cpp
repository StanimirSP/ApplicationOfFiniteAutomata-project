#include "constants.h"
#include <limits>
#include <stdexcept>

namespace Constants
{
	namespace Internal
	{
		static unsigned char SymbolTable[std::numeric_limits<unsigned char>::max()];
		constexpr unsigned char 	special_mask = 0b0000'0001;
		constexpr unsigned char	   operator_mask = 0b0000'0010;
		constexpr unsigned char parenthesis_mask = 0b0000'0100;
		constexpr unsigned char   forbidden_mask = 0b0000'1000;

		bool initTable()
		{
			// special symbols
			SymbolTable[static_cast<USymbol>(BaseElementBegin)] |= special_mask;
			SymbolTable[static_cast<USymbol>(BaseElementEnd)] |= special_mask;
			SymbolTable[static_cast<USymbol>(Union)] |= special_mask;
			SymbolTable[static_cast<USymbol>(Concatenation)] |= special_mask;
			SymbolTable[static_cast<USymbol>(KleeneStar)] |= special_mask;
			SymbolTable[static_cast<USymbol>(OpenParenthesis)] |= special_mask;
			SymbolTable[static_cast<USymbol>(CloseParenthesis)] |= special_mask;
			SymbolTable[static_cast<USymbol>(BasePlaceholder)] |= special_mask;
			SymbolTable[static_cast<USymbol>(EmptySet)] |= special_mask;

			// operators
			SymbolTable[static_cast<USymbol>(Union)] |= operator_mask;
			SymbolTable[static_cast<USymbol>(Concatenation)] |= operator_mask;
			SymbolTable[static_cast<USymbol>(KleeneStar)] |= operator_mask;

			// parentheses
			SymbolTable[static_cast<USymbol>(OpenParenthesis)] |= parenthesis_mask;
			SymbolTable[static_cast<USymbol>(CloseParenthesis)] |= parenthesis_mask;

			// symbols that should not appear in the alphabet
			SymbolTable[static_cast<USymbol>(ReplacementPos)] |= forbidden_mask;
			SymbolTable[static_cast<USymbol>(ReplacementStart)] |= forbidden_mask;
			SymbolTable[static_cast<USymbol>(ReplacementEnd)] |= forbidden_mask;
			SymbolTable[static_cast<USymbol>(Epsilon)] |= forbidden_mask;
			return true;
		}
		static const bool invokeInit = initTable();
	}
	bool isSpecial(USymbol c)
	{
		return Internal::SymbolTable[c] & Internal::special_mask;
	}
	bool isOperator(USymbol c)
	{
		return Internal::SymbolTable[c] & Internal::operator_mask;
	}
	bool isParenthesis(USymbol c)
	{
		return Internal::SymbolTable[c] & Internal::parenthesis_mask;
	}
	bool isForbidden(USymbol c)
	{
		return Internal::SymbolTable[c] & Internal::forbidden_mask;
	}
	int precedence(Symbol c)
	{
		switch(c)
		{
		case 		 Constants::Union: return 1;
		case Constants::Concatenation: return 2;
		case 	Constants::KleeneStar: return 3;
		default: throw std::logic_error("this should not happen");
		}
	}
}
