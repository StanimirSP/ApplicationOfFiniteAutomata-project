#ifndef MONOIDALFSA_H
#define MONOIDALFSA_H

#include <cstddef>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <queue>
#include <functional>
#include <utility>
#include <map>
#include <set>
#include <span>
#include <ranges>
#include <concepts>
#include <algorithm>
#include <array>
#include "constants.h"
#include "transition.h"
#include "regularExpression.h"

template<class LabelType>
class MonoidalFSA;

template<class LabelType>
MonoidalFSA<LabelType> regexToMFSA(RegularExpression<LabelType> re, const std::string& alphabet);

template<class LabelType>
class MonoidalFSA
{
protected:
	State statesCnt = 0;
	TransitionList<LabelType> transitions{statesCnt};
	std::unordered_set<State> initial, final;
	std::vector<Symbol> alphabet;
	std::array<int, std::numeric_limits<USymbol>::max() + 1> alphabetOrder;
	MonoidalFSA()
	{
		std::ranges::fill(alphabetOrder, -1);
	}
	friend MonoidalFSA<LabelType> regexToMFSA<>(RegularExpression<LabelType> re, const std::string& alphabet);
	friend class Bimachine;

	template<std::invocable<State> UnaryFunction, std::predicate<Transition<LabelType>> UnaryPredicate = Constants::AlwaysTrue>
	void BFS(State st, UnaryFunction use, UnaryPredicate pred = {}) const
	{
		std::queue<State> q;
		q.push(st);
		std::unordered_set<State> visited;
		visited.insert(st);
		while(!q.empty())
		{
			State currSt = q.front();
			q.pop();
			use(currSt);
			for(const auto& tr : transitions(currSt))
				if(pred(tr) && visited.insert(tr.To()).second)
					q.push(tr.To());
		}
	}
	void epsCloseInitial()
	{
		std::vector<State> epsClosure;
		auto appendToClosure = [&epsClosure](State st) { epsClosure.push_back(st); };
		auto isEpsTransition = [](const Transition<LabelType>& tr) { return tr.Label() == LabelType::epsilon(); };
		for(State init : initial)
			BFS(init, appendToClosure, isEpsTransition);
		initial.insert(epsClosure.begin(), epsClosure.end());
	}
	bool containsFinalState(const std::set<State>& set) const
	{
		for(State st : set)
			if(final.contains(st))
				return true;
		return false;
	}
	// precondition: transitions must be sorted by Label
	virtual std::vector<LabelType> findPseudoAlphabet() const
	{
		if(transitions.buffer.empty()) return {};
		std::vector<LabelType> pseudoAlphabet{transitions.buffer.front().Label()};
		for(const auto& tr : transitions.buffer)
			if(pseudoAlphabet.back() != tr.Label())
				pseudoAlphabet.push_back(tr.Label());
		return pseudoAlphabet;
	}
	// precondition: transitions must be sorted by From, then by Label
	MonoidalFSA& complete(const std::vector<LabelType>& pseudoAlphabet)
	{
		if(transitions.buffer.size() < statesCnt * pseudoAlphabet.size())
		{
			std::size_t trInd = 0, oldTrSize = transitions.buffer.size();
			for(State st = 0; st <= statesCnt; st++)
				for(const LabelType& letter : pseudoAlphabet)
					if(trInd >= oldTrSize || transitions.buffer[trInd].From() != st || transitions.buffer[trInd].Label() != letter)
						transitions.buffer.emplace_back(st, letter, statesCnt);
					else
						trInd++;
			statesCnt++;
		}
		return *this;
	}
	static std::unordered_set<State> filterAndRemap(const std::unordered_set<State>& states, const std::vector<State>& map)
	{
		std::unordered_set<State> remapped;
		for(State st : states)
			if(map[st] != Constants::InvalidState)
				remapped.insert(map[st]);
		return remapped;
	}
	void alphabetUnion(Symbol s)
	{
		if(alphabetOrder[static_cast<USymbol>(s)] < 0)
		{
			alphabetOrder[static_cast<USymbol>(s)] = alphabet.size();
			alphabet.push_back(s);
		}
	}
	void alphabetUnion(const MonoidalFSA& rhs)
	{
		for(Symbol s : rhs.alphabet)
			alphabetUnion(s);
	}
public:
	MonoidalFSA& removeEpsilon()
	{
		transitions.sort();
		epsCloseInitial();
		TransitionList<LabelType> newTransitions{statesCnt};
		newTransitions.startInd.reserve(statesCnt + 1);
		newTransitions.buffer.reserve(transitions.buffer.size());
		State lastState = Constants::InvalidState;
		auto isEpsTransition = [](const Transition<LabelType>& tr) { return tr.Label() == LabelType::epsilon(); };
		for(const auto& tr : transitions.buffer)
			if(tr.Label() != LabelType::epsilon())
			{
				newTransitions.startInd.insert(newTransitions.startInd.end(), tr.From() - lastState, newTransitions.buffer.size());
				lastState = tr.From();
				auto appendToNewTransitions = [&newTransitions, lastState, &lbl = tr.Label()](State dest) {
					newTransitions.buffer.emplace_back(lastState, lbl, dest);
					};
				BFS(tr.To(), appendToNewTransitions, isEpsTransition);
			}
		newTransitions.startInd.insert(newTransitions.startInd.end(), statesCnt - lastState, newTransitions.buffer.size());
		newTransitions.isSorted = true;
		transitions = std::move(newTransitions);
		return *this;
	}
	MonoidalFSA& trim()
	{
		transitions.sort();
		std::vector<State> newStates(statesCnt);
		auto markReachable = [&newStates](State st) { newStates[st]++; };
		auto notVisited = [&newStates](const Transition<LabelType>& tr) -> bool { return newStates[tr.To()] == 0; };
		for(State init : initial)
			if(newStates[init] == 0)
				BFS(init, markReachable, notVisited);
		transitions.reverse().sort();
		auto visitedOnce = [&newStates](const Transition<LabelType>& tr) -> bool { return newStates[tr.To()] == 1; };
		for(State fin : final)
			if(newStates[fin] == 1)
				BFS(fin, markReachable, visitedOnce);

		// rename states
		statesCnt = 0;
		for(State& st : newStates)
			st = (st == 2 ? statesCnt++ : Constants::InvalidState);

		// filter, remap and reverse transitions
		transitions.isSorted = false;
		std::size_t outputIndex = 0;
		for(auto& tr : transitions.buffer)
			if(newStates[tr.From()] != Constants::InvalidState && newStates[tr.To()] != Constants::InvalidState)
				transitions.buffer[outputIndex++] = {newStates[tr.To()], tr.Label(), newStates[tr.From()]};
		transitions.buffer.erase(transitions.buffer.begin() + outputIndex, transitions.buffer.end());

		//filter and remap initial and final states
		initial = filterAndRemap(initial, newStates);
		final = filterAndRemap(final, newStates);
		return *this;
	}
	template<std::invocable<LabelType, LabelType> Proj, std::predicate<LabelType, LabelType> Cond = Constants::AlwaysTrue>
	auto product(MonoidalFSA& rhs, Proj tranformLabel, Cond labelCondition = {})
	{
		this->transitions.sort();
		rhs.transitions.sort();
		MonoidalFSA<decltype(tranformLabel(std::declval<LabelType>(), std::declval<LabelType>()))> prod;

		// initialize alphabet
		prod.alphabetUnion(*this);
		prod.alphabetUnion(rhs);

		std::map<std::pair<State, State>, State> newStates;
		std::queue<std::pair<State, State>> q;
		// initialize newStates with I1 x I2
		for(State s1 : initial)
			for(State s2 : rhs.initial)
			{
				newStates[{s1, s2}] = prod.statesCnt++;
				q.emplace(s1, s2);
			}

		// initialize transitions
		prod.transitions.startInd.push_back(0);
		for(std::size_t step = 0; !q.empty(); step++)
		{
			auto [s1, s2] = q.front();
			q.pop();
			Transition<LabelType> eps1[1]{{s1, LabelType::epsilon(), s1}}, eps2[1]{{s2, LabelType::epsilon(), s2}};
			std::span<const Transition<LabelType>, std::dynamic_extent> span1{eps1}, span2{eps2};
			std::array<std::span<const Transition<LabelType>, std::dynamic_extent>, 2> r1{transitions(s1), span1}, r2{rhs.transitions(s2), span2};
			for(const auto& tr1 : r1 | std::views::join)
				for(const auto& tr2 : r2 | std::views::join)
					if(labelCondition(tr1.Label(), tr2.Label()))
					{
						auto [it, inserted] = newStates.try_emplace({tr1.To(), tr2.To()}, prod.statesCnt);
						if(inserted)
						{
							prod.statesCnt++;
							q.emplace(tr1.To(), tr2.To());
						}
						prod.transitions.buffer.emplace_back(step, tranformLabel(tr1.Label(), tr2.Label()), it->second); // it->second is the value for {tr1.To(), tr2.To()}
					}
			prod.transitions.startInd.push_back(prod.transitions.buffer.size());
		}
		prod.transitions.isSorted = true;

		// initialize I and F
		for(auto [st, newName] : newStates)
		{
			if(initial.contains(st.first) && rhs.initial.contains(st.second))
				prod.initial.insert(newName);
			if(final.contains(st.first) && rhs.final.contains(st.second))
				prod.final.insert(newName);
		}
		return prod.removeEpsilon().trim();
	}
	[[nodiscard]] MonoidalFSA Union(MonoidalFSA&& rhs) &&
	{
		if(transitions.buffer.size() < rhs.transitions.buffer.size())
			std::swap(*this, rhs);
		MonoidalFSA un(std::move(*this));
		un.transitions.isSorted = false;
		for(const auto& tr : rhs.transitions.buffer)
			un.transitions.buffer.emplace_back(tr.From() + un.statesCnt, tr.Label(), tr.To() + un.statesCnt);
		for(State init : rhs.initial)
			un.initial.insert(init + un.statesCnt);
		for(State fin : rhs.final)
			un.final.insert(fin + un.statesCnt);
		un.alphabetUnion(rhs);
		un.statesCnt += rhs.statesCnt;
		return un;
	}
	[[nodiscard]] MonoidalFSA Union(const MonoidalFSA& rhs) && { return std::move(*this).Union(MonoidalFSA{rhs}); }
	[[nodiscard]] MonoidalFSA Union(MonoidalFSA&& rhs) const& { return MonoidalFSA{*this}.Union(std::move(rhs)); }
	[[nodiscard]] MonoidalFSA Union(const MonoidalFSA& rhs) const& { return MonoidalFSA{*this}.Union(MonoidalFSA{rhs}); }
	[[nodiscard]] MonoidalFSA Concatenation(const MonoidalFSA& rhs) &&
	{
		MonoidalFSA concat;
		concat.transitions.buffer = std::move(transitions.buffer);
		for(const auto& tr : rhs.transitions.buffer)
			concat.transitions.buffer.emplace_back(tr.From() + statesCnt, tr.Label(), tr.To() + statesCnt);
		concat.initial = std::move(initial);
		for(State fin : rhs.final)
			concat.final.insert(fin + statesCnt);
		concat.alphabetUnion(*this);
		concat.alphabetUnion(rhs);
		for(State fin1 : final)
			for(State init2 : rhs.initial)
				concat.transitions.buffer.emplace_back(fin1, LabelType::epsilon(), init2 + statesCnt);
		concat.statesCnt = statesCnt + rhs.statesCnt;
		return concat;
	}
	[[nodiscard]] MonoidalFSA Concatenation(const MonoidalFSA& rhs) const& { return MonoidalFSA{*this}.Concatenation(rhs); }
	MonoidalFSA& Plus()
	{
		transitions.isSorted = false;
		for(State init : initial)
			transitions.buffer.emplace_back(statesCnt, LabelType::epsilon(), init);
		for(State fin : final)
			transitions.buffer.emplace_back(fin, LabelType::epsilon(), statesCnt);
		initial = {statesCnt++};
		return *this;
	}
	MonoidalFSA& KleeneStar()
	{
		Plus();
		final.insert(statesCnt - 1);
		return *this;
	}
	MonoidalFSA& Option()
	{
		initial.insert(statesCnt);
		final.insert(statesCnt++);
		transitions.startInd.push_back(transitions.buffer.size());
		return *this;
	}
	MonoidalFSA& pseudoDeterm()
	{
		this->removeEpsilon().trim().transitions.sort();
		std::map<std::set<State>, State> newStates;
		std::queue<const std::set<State>*> q;
		newStates[std::set<State>(this->initial.begin(), this->initial.end())] = 0;
		q.push(&newStates.begin()->first);
		this->statesCnt = 1;
		TransitionList<LabelType> newTransitions{this->statesCnt};
		std::unordered_set<State> newFinal;
		if(containsFinalState(*q.front()))
			newFinal.insert(0);
		std::map<LabelType, std::set<State>> nextSets;
		for(State step = 0; !q.empty(); step++)
		{
			auto& currSet = *q.front();
			q.pop();
			for(State st : currSet)
				for(const auto& tr : this->transitions(st))
					nextSets[tr.Label()].insert(tr.To());
			for(auto& label_set : nextSets)
			{
				auto [it, inserted] = newStates.try_emplace(std::move(label_set.second), this->statesCnt);
				if(inserted)
				{
					this->statesCnt++;
					q.push(&it->first);
					if(containsFinalState(it->first))
						newFinal.insert(it->second);
				}
				newTransitions.buffer.emplace_back(step, label_set.first, it->second);
			}
			nextSets.clear();
		}
		this->transitions = std::move(newTransitions);
		this->initial = {0};
		this->final = std::move(newFinal);
		return *this;
	}
	MonoidalFSA& pseudoMinimize()
	{
		pseudoDeterm();
		if(final.empty()) return *this;
		sortByLabel(transitions);
		std::vector<LabelType> pseudoAlphabet = findPseudoAlphabet();
		transitions.sort();
		complete(pseudoAlphabet);
		sortByLabel(transitions.reverse());
		transitions.sort();

		std::vector<State> classInd; // classInd[i] == j <=> state i belongs to equivalence class j
		classInd.reserve(statesCnt);

		auto cmpLabel = [](const Transition<LabelType>& a, const Transition<LabelType>& b) { return a.Label() < b.Label(); };
		struct EquivalenceClass
		{
			std::unordered_set<State> members; // which states are in the current class
			std::vector<State> plus; // subset of members
			std::vector<bool> inQueue; // inQueue[i] <=> (members, pseudoAlphabet[i]) in q
			EquivalenceClass(std::size_t alphabet_size): inQueue(alphabet_size) {}
		};
		std::vector<EquivalenceClass> eqClassBuffer;
		{ // initialize equivance classes
			EquivalenceClass finals(pseudoAlphabet.size()), nonfinals(pseudoAlphabet.size());
			for(State st = 0; st < statesCnt; st++)
				if(final.contains(st))
				{
					finals.members.insert(st);
					classInd.push_back(0);
				}
				else
				{
					nonfinals.members.insert(st);
					classInd.push_back(1);
				}
			eqClassBuffer.push_back(std::move(finals));
			if(!nonfinals.members.empty())
				eqClassBuffer.push_back(std::move(nonfinals));
		}

		std::queue<std::pair<State, State>> q; // queue of pairs ("equvance class index", "letter index")
		{ // initialize queue
			State smallerClassInd = (eqClassBuffer.size() < 2 || eqClassBuffer[0].members.size() < eqClassBuffer[1].members.size()) ? 0 : 1;
			for(std::size_t letterInd = 0; letterInd < pseudoAlphabet.size(); letterInd++)
			{
				q.emplace(smallerClassInd, letterInd);
				eqClassBuffer[smallerClassInd].inQueue[letterInd] = true;
			}
		}
		std::queue<State> splitQueue; // queue of class indices which may split; it is defined outside the outer loop to avoid reallocations
		while(!q.empty())
		{
			auto [classIndex, letterInd] = q.front();
			q.pop();
			eqClassBuffer[classIndex].inQueue[letterInd] = false;
			for(State s1 : eqClassBuffer[classIndex].members)
				for(const auto& tr : std::ranges::equal_range(transitions(s1), Transition{s1, pseudoAlphabet[letterInd], Constants::InvalidState}, cmpLabel))
				{ // (s1, a, tr.To()) is in reversed transitions
					State targetClassInd = classInd[tr.To()];
					EquivalenceClass& targetClass = eqClassBuffer[targetClassInd];
					if(targetClass.plus.empty())
						splitQueue.push(targetClassInd);
					targetClass.plus.push_back(tr.To());
				}

			while(!splitQueue.empty())
			{
				State toSplitInd = splitQueue.front();
				splitQueue.pop();
				if(eqClassBuffer[toSplitInd].members.size() != eqClassBuffer[toSplitInd].plus.size()) // toSplit indeed splits
				{
					eqClassBuffer.emplace_back(pseudoAlphabet.size());
					State newClassInd = eqClassBuffer.size() - 1;
					EquivalenceClass& newClass = eqClassBuffer[newClassInd];
					for(State toMove : eqClassBuffer[toSplitInd].plus)
					{
						eqClassBuffer[toSplitInd].members.erase(toMove);
						newClass.members.insert(toMove);
						classInd[toMove] = newClassInd;
					}
					for(std::size_t letterInd = 0; letterInd < pseudoAlphabet.size(); letterInd++)
						if(eqClassBuffer[toSplitInd].inQueue[letterInd])
						{
							q.emplace(newClassInd, letterInd);
							newClass.inQueue[letterInd] = true;
						}
						else
						{
							State smallerClassInd = eqClassBuffer[toSplitInd].members.size() < newClass.members.size() ? toSplitInd : newClassInd;
							q.emplace(smallerClassInd, letterInd);
							eqClassBuffer[smallerClassInd].inQueue[letterInd] = true;
						}
				}
				eqClassBuffer[toSplitInd].plus.clear();
			}
			// splitQueue is now empty and prepared for the next iteration
		}

		statesCnt = eqClassBuffer.size();
		initial = {classInd[*initial.begin()]}; // new initial state is the class of the old initial state
		{
			std::unordered_set<State> newFinal;
			for(State fin : final)
				newFinal.insert(classInd[fin]);
			final = std::move(newFinal);
		}
		TransitionList<LabelType> newTransitions{statesCnt};
		newTransitions.buffer.reserve(transitions.buffer.size());
		for(auto& tr : transitions.buffer) // transitions are still reversed!
			newTransitions.buffer.emplace_back(classInd[tr.To()], std::move(tr.Label()), classInd[tr.From()]);
		sortByLabel(newTransitions);
		newTransitions.sort();
		newTransitions.buffer.erase(std::ranges::unique(newTransitions.buffer).begin(), newTransitions.buffer.end());
		newTransitions.isSorted = false; // because newTransitions.startInd is now invalid
		transitions = std::move(newTransitions);
		return trim();
	}
	MonoidalFSA(const MonoidalFSA& rhs): statesCnt(rhs.statesCnt), transitions(rhs.transitions),
		initial(rhs.initial), final(rhs.final), alphabet(rhs.alphabet), alphabetOrder(rhs.alphabetOrder)
	{
		transitions.statesCnt = &statesCnt;
	}
	MonoidalFSA(MonoidalFSA&& rhs): statesCnt(std::exchange(rhs.statesCnt, 0)), transitions(std::move(rhs.transitions)),
		initial(std::move(rhs.initial)), final(std::move(rhs.final)), alphabet(std::move(rhs.alphabet)), alphabetOrder(std::move(rhs.alphabetOrder))
	{
		transitions.statesCnt = &statesCnt;
	}
	MonoidalFSA& operator=(MonoidalFSA rhs)
	{
		using std::swap;
		swap(statesCnt, rhs.statesCnt);
		swap(transitions, rhs.transitions);
		transitions.statesCnt = &statesCnt;
		swap(initial, rhs.initial);
		swap(final, rhs.final);
		swap(alphabet, rhs.alphabet);
		swap(alphabetOrder, rhs.alphabetOrder);
		return *this;
	}
	std::ostream& print(std::ostream& os = std::cout) const
	{
		os << "States: " << statesCnt << '\n';
		/*os << "Initial: ";
		for(State st : initial)
			os << st << ' ';
		os << "\nFinal: ";
		for(State st : final)
			os << st << ' ';*/
		os << "\nTransitions: " << transitions.buffer.size() << '\n';
		/*for(auto&& x : transitions.buffer)
			os << x << '\n';*/
		return os;
	}
	void clear()
	{
		statesCnt = 0;
		initial.clear();
		final.clear();
		transitions.clear();
		alphabet.clear();
		std::ranges::fill(alphabetOrder, -1);
	}
	friend std::ostream& operator<<(std::ostream& os, const MonoidalFSA& T)
	{
		for(Symbol s : T.alphabet)
			os << s;
		os << Constants::Epsilon << '\n'; // epsilon used as terminator
		os << T.statesCnt << ' ' << T.initial.size() << ' ' << T.final.size() << '\n';
		for(State st : T.initial)
			os << st << ' ';
		os << '\n';
		for(State st : T.final)
			os << st << ' ';
		return os << '\n' << T.transitions;
	}
	friend std::istream& operator>>(std::istream& is, MonoidalFSA& T)
	{
		T.clear();
		std::string alphabet;
		std::getline(is, alphabet, Constants::Epsilon);
		for(Symbol s : alphabet)
			T.alphabetUnion(s);

		std::size_t initialCnt, finalCnt;
		is >> T.statesCnt >> initialCnt >> finalCnt;
		while(initialCnt--)
		{
			State init;
			is >> init;
			T.initial.insert(init);
		}
		while(finalCnt--)
		{
			State fin;
			is >> fin;
			T.final.insert(fin);
		}
		return is >> T.transitions;
	}
};

#endif
