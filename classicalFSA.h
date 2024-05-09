#ifndef CLASSICALFSA_H
#define CLASSICALFSA_H

#include <concepts>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <utility>
#include <limits>
#include <cstddef>
#include <vector>
#include "monoidalFSA.h"
#include "transition.h"
#include "constants.h"

class ClassicalFSA: public MonoidalFSA<SymbolOrEpsilon>
{
	template<bool, class>
	friend class Transducer;
	friend class Bimachine;
	ClassicalFSA() = default;

	virtual std::vector<SymbolOrEpsilon> findPseudoAlphabet() const override
	{
		return {alphabet.begin(), alphabet.end()};
	}
	std::map<std::set<State>, State> convertToDFSA_ret(bool noTrim = false)
	{
		if(!noTrim)
			this->removeEpsilon().trim();
		this->transitions.sort();
		std::map<std::set<State>, State> newStates;
		std::queue<const std::set<State>*> q;
		newStates[std::set<State>(this->initial.begin(), this->initial.end())] = 0;
		q.push(&newStates.begin()->first);
		this->statesCnt = 1;
		TransitionList<SymbolOrEpsilon> newTransitions{this->statesCnt};
		std::unordered_set<State> newFinal;
		if(containsFinalState(*q.front()))
			newFinal.insert(0);
		std::set<State> nextSets[std::numeric_limits<USymbol>::max() + 1];
		newTransitions.startInd.push_back(0);
		for(State step = 0; !q.empty(); step++)
		{
			auto& currSet = *q.front();
			q.pop();
			for(State st : currSet)
				for(auto& tr : this->transitions(st))
					nextSets[static_cast<USymbol>(tr.Label())].insert(tr.To());
			for(USymbol letter : this->alphabet)
			{
				auto [it, inserted] = newStates.try_emplace(std::move(nextSets[letter]), this->statesCnt);
				nextSets[letter].clear();
				if(inserted)
				{
					this->statesCnt++;
					q.push(&it->first);
					if(containsFinalState(it->first))
						newFinal.insert(it->second);
				}
				newTransitions.buffer.emplace_back(step, letter, it->second);
			}
			newTransitions.startInd.push_back(newTransitions.buffer.size());
		}
		newTransitions.isSorted = true;
		this->transitions = std::move(newTransitions);
		this->initial = {0};
		this->final = std::move(newFinal);
		return newStates;
	}
public:
	ClassicalFSA(const MonoidalFSA<SymbolOrEpsilon>& mFSA): MonoidalFSA<SymbolOrEpsilon>(mFSA) {}
	ClassicalFSA(MonoidalFSA<SymbolOrEpsilon>&& mFSA): MonoidalFSA<SymbolOrEpsilon>(std::move(mFSA)) {}

	ClassicalFSA& convertToDFSA()
	{
		convertToDFSA_ret();
		return *this;
	}
	ClassicalFSA& complement()
	{
		convertToDFSA();
		std::unordered_set<State> newFinal;
		for(State i = 0; i < this->statesCnt; i++)
			if(!this->final.contains(i))
				newFinal.insert(i);
		this->final = std::move(newFinal);
		return *this;
	}
	ClassicalFSA intersect(ClassicalFSA& rhs)
	{
		convertToDFSA();
		rhs.convertToDFSA();
		auto tranformLabel = [](SymbolOrEpsilon first, SymbolOrEpsilon second) -> SymbolOrEpsilon {
			return first;
		};
		auto labelCondition = [](SymbolOrEpsilon first, SymbolOrEpsilon second) -> bool {
			return first == second && first != Constants::Epsilon;
		};
		return this->product(rhs, tranformLabel, labelCondition);
	}
	static ClassicalFSA createFromSymbolSet(const std::ranges::input_range auto& symbols)
	{
		ClassicalFSA res;
		res.initial.insert(res.statesCnt++);
		for(Symbol c : symbols)
		{
			res.alphabetUnion(c);
			res.transitions.buffer.emplace_back(0, c, res.statesCnt);
			res.final.insert(res.statesCnt++);
		}
		return res;
	}
};

#endif
