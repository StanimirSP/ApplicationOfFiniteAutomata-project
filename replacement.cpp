#include <concepts>
#include <stdexcept>
#include <string>
#include <vector>
#include "transducer.h"
#include "classicalFSA.h"
#include "constants.h"
#include "regularExpression.h"
#include "replacement.h"

namespace Internal
{
	std::string AlphaCompl(const std::string& alpha, const std::string& extendedAlphabet)
	{
		std::string alphaCompl;
		for(Symbol s : extendedAlphabet)
			if(!alpha.contains(s))
				alphaCompl.push_back(s);
		return alphaCompl;
	}
	LetterTransducer intro(const std::string& what, const std::string& extendedAlphabet)
	{
		std::vector<SymbolPair> eps_what;
		eps_what.reserve(what.size());
		for(Symbol s : what)
			eps_what.push_back({Constants::Epsilon, s});
		LetterTransducer id = LetterTransducer::identity(ClassicalFSA::createFromSymbolSet(AlphaCompl(what, extendedAlphabet)));
		LetterTransducer intr = LetterTransducer::createFromSymbolSet(eps_what);
		return std::move(id).Union(std::move(intr)).KleeneStar();
	}
	LetterTransducer introX(const std::string& what, const std::string& extendedAlphabet)
	{
		LetterTransducer id = LetterTransducer::identity(ClassicalFSA::createFromSymbolSet(AlphaCompl(what, extendedAlphabet)));
		return intro(what, extendedAlphabet).Concatenation(std::move(id)).Option();
	}
	LetterTransducer Xintro(const std::string& what, const std::string& extendedAlphabet)
	{
		LetterTransducer id = LetterTransducer::identity(ClassicalFSA::createFromSymbolSet(AlphaCompl(what, extendedAlphabet)));
		return std::move(id).Concatenation(intro(what, extendedAlphabet)).Option();
	}
	ClassicalFSA ignX(LetterTransducer& identityL, const std::string& what, const std::string& extendedAlphabet)
	{
		LetterTransducer introx = introX(what, extendedAlphabet).pseudoMinimize();
		return identityL.compose(introx).Range();
	}
	ClassicalFSA Xign(LetterTransducer& identityL, const std::string& what, const std::string& extendedAlphabet)
	{
		LetterTransducer xintro = Xintro(what, extendedAlphabet).pseudoMinimize();
		return identityL.compose(xintro).Range();
	}
	ClassicalFSA PiffS(ClassicalFSA&& P, ClassicalFSA&& S, const std::string& extendedAlphabet)
	{
		ClassicalFSA ifPthenS = P.Concatenation(ClassicalFSA(S).complement());
		ClassicalFSA ifSthenP = std::move(P.complement()).Concatenation(S);
		return ifPthenS.complement().intersect(ifSthenP.complement());
	}
	ClassicalFSA LiffR(ClassicalFSA&& L, ClassicalFSA&& R, const std::string& extendedAlphabet)
	{
		ClassicalFSA ALL = ClassicalFSA::createFromSymbolSet(extendedAlphabet).KleeneStar();
		return PiffS(ALL.Concatenation(std::move(L)), R.Concatenation(ALL), extendedAlphabet);
	}
	LetterTransducer InitialMatch(LetterTransducer& IdDomT, const std::string& extendedAlphabet)
	{
		const std::string replPos{Constants::ReplacementPos};
		LetterTransducer lhs = intro(replPos, extendedAlphabet).pseudoMinimize();
		LetterTransducer rhs = LetterTransducer::identity(LiffR(ClassicalFSA::createFromSymbolSet(replPos), Xign(IdDomT, replPos, extendedAlphabet), extendedAlphabet).pseudoMinimize());
		return lhs.compose(rhs);
	}
	LetterTransducer LongestMatch(LetterTransducer& IdDomT, const std::string& extendedAlphabet)
	{
		ClassicalFSA containsReplEnd = contain(ClassicalFSA::createFromSymbolSet(std::string{Constants::ReplacementEnd}), extendedAlphabet).pseudoMinimize();
		return LetterTransducer::identity(contain(
			ClassicalFSA::createFromSymbolSet(std::string{Constants::ReplacementStart})
			.Concatenation(ignX(IdDomT, {Constants::ReplacementStart, Constants::ReplacementEnd}, extendedAlphabet).intersect(containsReplEnd)).pseudoMinimize(),
			extendedAlphabet
		).complement());
	}
}
