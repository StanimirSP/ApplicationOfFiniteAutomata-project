#ifndef REPLACEMENT_H
#define REPLACEMENT_H

#include <concepts>
#include <stdexcept>
#include <string>
#include <vector>
#include "transducer.h"
#include "classicalFSA.h"
#include "constants.h"
#include "regularExpression.h"

LetterTransducer Replace(LetterTransducer&& T, const std::ranges::input_range auto& alphabet);
LetterTransducer Replace(const LetterTransducer& T, const std::ranges::input_range auto& alphabet);

namespace Internal
{
	ClassicalFSA contain(ClassicalFSA&& L, const std::ranges::input_range auto& alphabet)
	{
		ClassicalFSA all = ClassicalFSA::createFromSymbolSet(alphabet).KleeneStar();
		return all.Concatenation(std::move(L)).Concatenation(std::move(all));
	}
	LetterTransducer NT(ClassicalFSA&& domainT, const std::ranges::input_range auto& alphabet)
	{
		ClassicalFSA tmp = contain(std::move(domainT), alphabet);
		tmp.complement().trim();
		return LetterTransducer::identity(std::move(tmp));
	}
	std::string createExtendedAlphabet(const std::ranges::input_range auto& alphabet)
	{
		std::string extendedAlphabet;
		if constexpr(std::ranges::sized_range<decltype(alphabet)>)
			extendedAlphabet.reserve(alphabet.size() + 3);
		for(Symbol s : alphabet)
			if(!Constants::isForbidden(s)) extendedAlphabet.push_back(s);
			else throw std::runtime_error("alphabet contains forbidden symbols");
		extendedAlphabet.push_back(Constants::ReplacementPos);
		extendedAlphabet.push_back(Constants::ReplacementStart);
		extendedAlphabet.push_back(Constants::ReplacementEnd);
		return extendedAlphabet;
	}
	std::string AlphaCompl(const std::string& alpha, const std::string& extendedAlphabet);
	LetterTransducer intro(const std::string& what, const std::string& extendedAlphabet);
	LetterTransducer introX(const std::string& what, const std::string& extendedAlphabet);
	LetterTransducer Xintro(const std::string& what, const std::string& extendedAlphabet);
	ClassicalFSA ignX(LetterTransducer& identityL, const std::string& what, const std::string& extendedAlphabet);
	ClassicalFSA Xign(LetterTransducer& identityL, const std::string& what, const std::string& extendedAlphabet);
	ClassicalFSA PiffS(ClassicalFSA&& P, ClassicalFSA&& S, const std::string& extendedAlphabet);
	ClassicalFSA LiffR(ClassicalFSA&& L, ClassicalFSA&& R, const std::string& extendedAlphabet);
	LetterTransducer InitialMatch(LetterTransducer& IdDomT, const std::string& extendedAlphabet);
	LetterTransducer LeftToRight(LetterTransducer& IdDomT, const std::ranges::input_range auto& alphabet, const std::string& extendedAlphabet)
	{
		const std::string replPos{Constants::ReplacementPos};
		LetterTransducer idAll = LetterTransducer::identity(ClassicalFSA::createFromSymbolSet(alphabet).KleeneStar());
		LetterTransducer ltr = idAll
			.Concatenation(LetterTransducer::createFromSymbolSet(std::vector{SymbolPair{Constants::ReplacementPos, Constants::ReplacementStart}}))
			.Concatenation(LetterTransducer::identity(ignX(IdDomT, replPos, extendedAlphabet)))
			.Concatenation(LetterTransducer::createFromSymbolSet(std::vector{SymbolPair{Constants::Epsilon, Constants::ReplacementEnd}}))
			.KleeneStar().Concatenation(std::move(idAll)).pseudoMinimize();
		LetterTransducer removeReplPos = Replace(LetterTransducer::createFromSymbolSet(std::vector{SymbolPair{Constants::ReplacementPos, Constants::Epsilon}}), extendedAlphabet).pseudoMinimize();
		return ltr.compose(removeReplPos);
	}
	LetterTransducer LongestMatch(LetterTransducer& IdDomT, const std::string& extendedAlphabet);
}

LetterTransducer Replace(LetterTransducer&& T, const std::ranges::input_range auto& alphabet)
{
	auto NT = Internal::NT(T.Domain(), alphabet);
	return NT.Concatenation(std::move(T)).KleeneStar().Concatenation(std::move(NT));
}

LetterTransducer Replace(const LetterTransducer& T, const std::ranges::input_range auto& alphabet)
{
	return Replace(LetterTransducer{T}, alphabet);
}

LetterTransducer ReplaceLML(LetterTransducer&& T, const std::ranges::input_range auto& alphabet)
{
	if(T.epsilonInDom())
		throw std::runtime_error("Dom(T) contains epsilon");
	std::string extendedAlphabet = Internal::createExtendedAlphabet(alphabet);
	LetterTransducer IdDomT = LetterTransducer::identity(T.Domain()).pseudoMinimize();
	LetterTransducer initialMatch = Internal::InitialMatch(IdDomT, extendedAlphabet).pseudoMinimize();
	LetterTransducer leftToRight = Internal::LeftToRight(IdDomT, alphabet, extendedAlphabet).pseudoMinimize();
	LetterTransducer longestMatch = Internal::LongestMatch(IdDomT, extendedAlphabet).pseudoMinimize();
	LetterTransducer replacement = Replace(
		LetterTransducer::createFromSymbolSet(std::vector{SymbolPair{Constants::ReplacementStart, Constants::Epsilon}})
			.Concatenation(std::move(T))
			.Concatenation(LetterTransducer::createFromSymbolSet(std::vector{SymbolPair{Constants::ReplacementEnd, Constants::Epsilon}})),
		extendedAlphabet
	).pseudoMinimize();
	return initialMatch.compose(leftToRight).compose(longestMatch).compose(replacement);
}

LetterTransducer ReplaceLML(const LetterTransducer& T, const std::ranges::input_range auto& alphabet)
{
	return ReplaceLML(LetterTransducer{T}, alphabet);
}

#endif
