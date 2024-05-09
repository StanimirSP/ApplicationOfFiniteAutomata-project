#ifndef BIMACHINE_H
#define BIMACHINE_H

#include "classicalFSA.h"
#include "transducer.h"
#include "transition.h"
#include "constants.h"
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <queue>
#include <utility>
#include <iostream>
#include <iomanip>

struct BimachineOutputElement
{
	State left;
	Symbol label;
	State right;
	Word out;
	constexpr BimachineOutputElement() = default;
	BimachineOutputElement(State left, Symbol label, State right, const Word& out): left(left), label(label), right(right)
	{
		for(Symbol s : out)
			if(s != Constants::Epsilon)
				this->out.push_back(s);
	}
	friend std::ostream& operator<<(std::ostream& os, const BimachineOutputElement& bmo)
	{
		return os << bmo.left << ' ' << bmo.label << ' ' << bmo.right << ' ' << std::quoted(bmo.out);
	}
	friend std::istream& operator>>(std::istream& is, BimachineOutputElement& bmo)
	{
		return is >> bmo.left >> bmo.label >> bmo.right >> std::quoted(bmo.out);
	}
};

struct BimachineOutput
{
	friend class Bimachine;
	std::vector<BimachineOutputElement> buffer;
	const State* rightStatesCnt;

	BimachineOutput(const State& rightStatesCnt):
		rightStatesCnt(&rightStatesCnt) {}
	void sort()
	{
		// already sorted by left by the construction
		sort(std::numeric_limits<std::make_unsigned_t<Symbol>>::max(), &BimachineOutputElement::label);
		sort(static_cast<std::size_t>(*rightStatesCnt) - 1, &BimachineOutputElement::right);
	}
	const Word& operator()(State left, Symbol label, State right) const
	{
		const BimachineOutputElement value(left, label, right, Word{});
		auto cmpLabel = [](const BimachineOutputElement& a, const BimachineOutputElement& b) { return a.label < b.label; };
		auto Right_Label = std::ranges::equal_range((*this)(right), value, cmpLabel);
		auto cmpLeft = [](const BimachineOutputElement& a, const BimachineOutputElement& b) { return a.left < b.left; };
		return std::ranges::lower_bound(Right_Label, value, cmpLeft)->out;
		/*auto Right_Label_Left = std::ranges::equal_range(Right_Label, value, cmpLeft);
		if(Right_Label_Left.size() < 1)
			throw std::runtime_error("very bad: output not defined but it should be");
		if(Right_Label_Left.size() > 1)
			throw std::runtime_error("very bad: output not unique but it should be");
		return Right_Label_Left.begin()->out;*/
	}
	void clear()
	{
		startInd.clear();
		buffer.clear();
	}
	friend std::ostream& operator<<(std::ostream& os, const BimachineOutput& bmo)
	{
		os << bmo.buffer.size() << '\n';
		for(auto&& el : bmo.buffer)
			os << el << '\n';
		os << bmo.startInd.size() << '\n';
		for(std::size_t ind : bmo.startInd)
			os << ind << ' ';
		return os << '\n';
	}
	friend std::istream& operator>>(std::istream& is, BimachineOutput& bmo)
	{
		bmo.clear();
		std::size_t functionSize;
		is >> functionSize;
		bmo.buffer.reserve(functionSize);
		while(functionSize--)
		{
			BimachineOutputElement el;
			is >> el;
			bmo.buffer.push_back(std::move(el));
		}
		std::size_t startIndSize;
		is >> startIndSize;
		bmo.startInd.reserve(startIndSize);
		while(startIndSize--)
		{
			std::size_t ind;
			is >> ind;
			bmo.startInd.push_back(ind);
		}
		return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
private:
	std::vector<std::size_t> startInd;
	template<std::invocable<BimachineOutputElement> Proj>
	void sort(std::size_t maxValue, Proj proj)
	{
		count(maxValue, proj);
		for(std::size_t i = 1; i < startInd.size(); i++)
			startInd[i] += startInd[i - 1];
		std::vector<BimachineOutputElement> sorted(buffer.size());
		for(auto& tr : buffer | std::views::reverse)
			sorted[--startInd[std::invoke(proj, tr)]] = std::move(tr);
		buffer = std::move(sorted);
	}
	void count(std::size_t maxValue, std::invocable<BimachineOutputElement> auto proj)
	{
		startInd.resize(maxValue + 2);
		std::ranges::fill(startInd, 0);
		for(const auto& tr : buffer)
			startInd[std::invoke(proj, tr)]++;
	}
	std::span<const BimachineOutputElement> operator()(State rightState) const
	{
		return {buffer.begin() + startInd[rightState], buffer.begin() + startInd[rightState + 1]};
	}
};

class Bimachine
{
	ClassicalFSA left, right;
	BimachineOutput output;

	struct LeftState
	{
		std::set<State> subset;
		std::map<State, State> function;
		void clear()
		{
			subset.clear();
			function.clear();
		}
		auto operator<=>(const LeftState&) const = default;
	};

	std::vector<std::set<State>> initializeRight(const Transducer<false, Symbol_Word>& rtime, std::string_view alphabet)
	{
		right.statesCnt = rtime.statesCnt;
		for(Symbol s : alphabet)
			right.alphabetUnion(s);
		right.initial = rtime.final; // for inverse determinization
		for(const auto& tr : rtime.transitions.buffer)
			right.transitions.buffer.emplace_back(tr.To(), tr.Label().first, tr.From()); // reverse transitions and project on tape 1
		std::map<std::set<State>, State> states = right.convertToDFSA_ret(true);
		std::vector<std::set<State>> ret(states.size());
		while(!states.empty())
		{
			auto handle = states.extract(states.begin());
			ret[handle.mapped()] = std::move(handle.key());
		}
		return ret;
	}
	LeftState initialLeft(const Transducer<false, Symbol_Word>& rtime, const std::vector<std::set<State>>& statesRight)
	{
		LeftState initialLeft;
		initialLeft.subset.insert(rtime.initial.begin(), rtime.initial.end());
		std::vector<State> intersection;
		for(State st = 0; st < statesRight.size(); st++)
		{
			std::ranges::set_intersection(statesRight[st], initialLeft.subset, std::back_inserter(intersection));
			if(!intersection.empty())
				initialLeft.function[st] = intersection.front(); // arbitrary element of the intersection
			intersection.clear();
		}
		return initialLeft;
	}
	std::vector<State> findRightPath(const Word& input) const
	{
		std::vector<State> path;
		path.reserve(input.size() + 1);
		State currSt = *right.initial.begin();
		path.push_back(currSt);
		for(Symbol s : input | std::views::reverse)
		{
			int letterInd = right.alphabetOrder[static_cast<USymbol>(s)];
			if(letterInd < 0)
				throw std::runtime_error("input string contains '" + std::string{s} + "' which is not in the alphabet");
			currSt = (right.transitions(currSt).begin() + letterInd)->To();
			path.push_back(currSt);
		}
		return path;
	}
public:
	Bimachine(const Transducer<false, Symbol_Word>& rtime, std::string_view alphabet): Bimachine(Transducer{rtime}, alphabet) {}
	Bimachine(Transducer<false, Symbol_Word>&& rtime, std::string_view alphabet): output(right.statesCnt)
	{
		std::vector<std::set<State>> statesRight = initializeRight(rtime, alphabet);
		std::map<LeftState, State> leftStates;
		leftStates[initialLeft(rtime, statesRight)] = left.statesCnt++;
		std::queue<const LeftState*> q;
		q.push(&leftStates.begin()->first);

		TransitionList<SymbolOrEpsilon> reversedRight = right.transitions;
		reversedRight.reverse().sort(); // these are actually in the original direction since the right automaton is reversed
		sortByLabel(rtime.transitions);
		rtime.transitions.sort();
		LeftState nextLeftStates[std::numeric_limits<USymbol>::max() + 1];
		auto cmpLabelFirst = [](const Transition<Symbol_Word>& a, const Transition<Symbol_Word>& b) { return a.Label().first < b.Label().first; };
		left.transitions.startInd.push_back(0);
		for(State step = 0; !q.empty(); step++)
		{
			const LeftState& currLeftSt = *q.front();
			q.pop();
			for(State st : currLeftSt.subset)
				for(const auto& tr : rtime.transitions(st))
					nextLeftStates[static_cast<USymbol>(tr.Label().first)].subset.insert(tr.To());
			for(auto [R, st] : currLeftSt.function)
				for(const auto& tr : reversedRight(R))
					for(const auto& rtimeTr :
						std::ranges::equal_range(rtime.transitions(st),
							Transition{Constants::InvalidState, Symbol_Word{tr.Label(), Word{}}, Constants::InvalidState}, cmpLabelFirst))
						if(statesRight[tr.To()].contains(rtimeTr.To()))
						{
							nextLeftStates[static_cast<USymbol>(tr.Label())].function[tr.To()] = rtimeTr.To();
							output.buffer.emplace_back(step, tr.Label(), tr.To(), rtimeTr.Label().second);
							break;
						}

			for(USymbol letter : right.alphabet)
			{
				auto [it, inserted] = leftStates.try_emplace(std::move(nextLeftStates[letter]), left.statesCnt);
				nextLeftStates[letter].clear();
				if(inserted)
				{
					left.statesCnt++;
					q.push(&it->first);
				}
				left.transitions.buffer.emplace_back(step, letter, it->second);
			}
			left.transitions.startInd.push_back(left.transitions.buffer.size());
		}
		left.transitions.isSorted = true;
		left.initial = {0};
		output.sort();

		// debug:
		/*std::cerr << left << '\n';
		std::cerr << right << '\n';
		std::cerr << "Output function size: " << output.buffer.size() << '\n';*/
		/*for(auto& [left, label, right, out] : output.buffer)
			std::cerr << left << ' ' << label << ' ' << right << " -> " << std::quoted(out) << '\n';*/
	}
	Bimachine(const Bimachine& rhs): left(rhs.left), right(rhs.right), output(rhs.output)
	{
		output.rightStatesCnt = &right.statesCnt;
	}
	Bimachine(Bimachine&& rhs): left(std::move(rhs.left)), right(std::move(rhs.right)), output(std::move(rhs.output))
	{
		output.rightStatesCnt = &right.statesCnt;
	}
	Bimachine& operator=(Bimachine rhs)
	{
		using std::swap;
		swap(left, rhs.left);
		swap(right, rhs.right);
		swap(output, rhs.output);
		output.rightStatesCnt = &right.statesCnt;
		return *this;
	}
	Word operator()(const Word& input) const
	{
		std::vector<State> path = findRightPath(input);
		Word out;
		State currLeftSt = *left.initial.begin();
		for(auto revPathIt = path.rbegin(); Symbol s : input)
		{
			out += output(currLeftSt, s, *++revPathIt);
			int letterInd = right.alphabetOrder[static_cast<USymbol>(s)]; // right alphabet because left alphabet is the same and there is no need to be stored
			currLeftSt = (left.transitions(currLeftSt).begin() + letterInd)->To();
		}
		return out;
	}
	friend std::ostream& operator<<(std::ostream& os, const Bimachine& bm)
	{
		return os << bm.left << bm.right << bm.output;
	}
	friend std::istream& operator>>(std::istream& is, Bimachine& bm)
	{
		return is >> bm.left >> bm.right >> bm.output;
	}
};

#endif
