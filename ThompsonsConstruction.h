#ifndef THOMPSONSCONSTRUCTION_H
#define THOMPSONSCONSTRUCTION_H

#include "constants.h"
#include "monoidalFSA.h"
#include "regularExpression.h"
#include <vector>
#include <stack>
#include <utility>

template<class LabelType>
struct ThompsonState;

template<class LabelType>
struct ThompsonTransition
{
	LabelType label;
	State indNext = Constants::InvalidState;
};

template<class LabelType>
struct ThompsonState
{
	ThompsonTransition<LabelType> transitions[2];
};

template<class LabelType>
class ThompsonAutomaton
{
	State indStart, indFinal;
	ThompsonAutomaton(State indStart, State indFinal): indStart(indStart), indFinal(indFinal) {}
public:
	static ThompsonAutomaton createForEmpty(std::vector<ThompsonState<LabelType>>& stateBuffer)
	{
		stateBuffer.resize(stateBuffer.size() + 2); // 2 new states
		return {static_cast<State>(stateBuffer.size() - 2), static_cast<State>(stateBuffer.size() - 1)};
	}
	static ThompsonAutomaton createForBase(LabelType&& label, std::vector<ThompsonState<LabelType>>& stateBuffer)
	{
		stateBuffer.emplace_back(); // new final state
		stateBuffer.push_back({{ {std::move(label), static_cast<State>(stateBuffer.size() - 1)} }}); // new initial state
		return {static_cast<State>(stateBuffer.size() - 1), static_cast<State>(stateBuffer.size() - 2)};
	}
	ThompsonAutomaton Union(const ThompsonAutomaton& rhs, std::vector<ThompsonState<LabelType>>& stateBuffer) const
	{
		stateBuffer.push_back({{ {LabelType::epsilon(), indStart}, {LabelType::epsilon(), rhs.indStart} }}); // new initial state
		stateBuffer[indFinal].	  transitions[0] = {LabelType::epsilon(), static_cast<State>(stateBuffer.size())};
		stateBuffer[rhs.indFinal].transitions[0] = {LabelType::epsilon(), static_cast<State>(stateBuffer.size())};
		stateBuffer.emplace_back(); // new final state
		return {static_cast<State>(stateBuffer.size() - 2), static_cast<State>(stateBuffer.size() - 1)};
	}
	ThompsonAutomaton Concatenation(const ThompsonAutomaton& rhs, std::vector<ThompsonState<LabelType>>& stateBuffer) const
	{
		stateBuffer[indFinal].transitions[0] = {LabelType::epsilon(), rhs.indStart};
		return {indStart, rhs.indFinal};
	}
	ThompsonAutomaton KleeneStar(std::vector<ThompsonState<LabelType>>& stateBuffer) const
	{
		stateBuffer.emplace_back(); // new final state
		stateBuffer[indFinal].transitions[0] = {LabelType::epsilon(), indStart};
		stateBuffer[indFinal].transitions[1] = {LabelType::epsilon(), static_cast<State>(stateBuffer.size() - 1)};
		stateBuffer.push_back({{ {LabelType::epsilon(), indStart}, {LabelType::epsilon(), static_cast<State>(stateBuffer.size() - 1)} }}); // new initial state
		return {static_cast<State>(stateBuffer.size() - 1), static_cast<State>(stateBuffer.size() - 2)};
	}
	std::vector<Transition<LabelType>> Transitions(std::vector<ThompsonState<LabelType>>& stateBuffer) const
	{
		std::vector<Transition<LabelType>> res;
		for(State i = 0; i < stateBuffer.size(); i++)
			for(auto& tr : stateBuffer[i].transitions)
				if(tr.indNext != Constants::InvalidState)
					res.emplace_back(i, std::move(tr.label), tr.indNext);
		return res;
	}
	friend MonoidalFSA<LabelType> regexToMFSA<>(RegularExpression<LabelType> re, const std::string& alphabet);
};

template<class LabelType>
ThompsonAutomaton<LabelType> ThompsonConstruction(RegularExpression<LabelType>&& re, std::vector<ThompsonState<LabelType>>& stateBuffer)
{
	stateBuffer.clear();
	stateBuffer.reserve(2 * re.TokenizedReversePolishNotation().size());
	std::stack<ThompsonAutomaton<LabelType>> st;
	auto baseElemIt = const_cast<std::vector<LabelType>&>(re.BaseTokens()).begin();
	for(Symbol c : re.TokenizedReversePolishNotation())
		switch(c)
		{
		case Constants::EmptySet:
			st.push(ThompsonAutomaton<LabelType>::createForEmpty(stateBuffer));
			break;
		case Constants::BasePlaceholder:
			st.push(ThompsonAutomaton<LabelType>::createForBase(std::move(*baseElemIt++), stateBuffer));
			break;
		case Constants::KleeneStar:
		{
			if(st.empty()) throw std::runtime_error("Bad regular expression");
			auto tmp = st.top().KleeneStar(stateBuffer);
			st.pop();
			st.push(tmp);
			break;
		}
		{
			ThompsonAutomaton<LabelType>(ThompsonAutomaton<LabelType>:: * concatOrUnoin)(const ThompsonAutomaton<LabelType>&, std::vector<ThompsonState<LabelType>>&) const;
		case Constants::Concatenation:
			concatOrUnoin = &ThompsonAutomaton<LabelType>::Concatenation;
			goto concat;
		case Constants::Union:
			concatOrUnoin = &ThompsonAutomaton<LabelType>::Union;
		concat:
			if(st.size() < 2) throw std::runtime_error("Bad regular expression");
			auto rhs = st.top();
			st.pop();
			auto lhs = st.top();
			st.pop();
			st.push((lhs.*concatOrUnoin)(rhs, stateBuffer));
			break;
		}
		default:
			throw std::runtime_error("Bad regular expression: unexpected character in tokenized RPN");
		}
	if(st.size() != 1) throw std::runtime_error("Bad regular expression");
	return st.top();
}

template<class LabelType>
MonoidalFSA<LabelType> regexToMFSA(RegularExpression<LabelType> re, const std::string& alphabet)
{
	std::vector<ThompsonState<LabelType>> stateBuffer;
	const ThompsonAutomaton<LabelType> Thompson = ThompsonConstruction(std::move(re), stateBuffer);
	MonoidalFSA<LabelType> res;
	res.statesCnt = stateBuffer.size();
	res.initial.insert(Thompson.indStart);
	res.final.insert(Thompson.indFinal);
	res.transitions.buffer = Thompson.Transitions(stateBuffer);
	for(Symbol c : alphabet)
		res.alphabetUnion(c);
	return res;
}

#endif
