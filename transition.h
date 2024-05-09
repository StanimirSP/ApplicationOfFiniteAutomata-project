#ifndef TRANSITION_H
#define TRANSITION_H

#include "constants.h"
#include "regularExpression.h"
#include <utility>
#include <vector>
#include <algorithm>
#include <span>
#include <iostream>
#include <concepts>
#include <functional>
#include <type_traits>
#include <limits>
#include <ranges>
#include <compare>

template<class LabelType>
class Transition
{
	State from;
	LabelType label;
	State to;
public:
	constexpr Transition() = default;
	Transition(State from, const LabelType& label, State to) noexcept: from(from), label(label), to(to) {}
	Transition(State from, LabelType&& label, State to) noexcept: from(from), label(std::move(label)), to(to) {}
	State From() const noexcept
	{
		return from;
	}
	const LabelType& Label() const noexcept
	{
		return label;
	}
	State To() const noexcept
	{
		return to;
	}
	void reverse() noexcept
	{
		std::swap(from, to);
	}
	auto operator<=>(const Transition& rhs) const = default;
	bool operator==(const Transition&) const = default;
	friend std::ostream& operator<<(std::ostream& os, const Transition& tr)
	{
		return os << tr.From() << ' ' << tr.Label() << ' ' << tr.To();
	}
	friend std::istream& operator>>(std::istream& is, Transition& tr)
	{
		return is >> tr.from >> tr.label >> tr.to;
	}
};

template<class LabelType>
struct TransitionList;

template<class LabelType>
void sortByLabel(TransitionList<LabelType>& list);

template<class LabelType>
struct TransitionList
{
	const State* statesCnt;
	std::vector<std::size_t> startInd;
	std::vector<Transition<LabelType>> buffer;
	bool isSorted = false; // if false, startInd is invalid

	TransitionList(const State& statesCnt): statesCnt(&statesCnt) {}
	void sort()
	{
		if(!isSorted) sort(static_cast<std::size_t>(*statesCnt) - 1, &Transition<LabelType>::From);
		isSorted = true;
	}
	TransitionList& reverse() noexcept
	{
		isSorted = false;
		for(auto& tr : buffer)
			tr.reverse();
		return *this;
	}
	TransitionList& reflexiveClose()
	{
		isSorted = false;
		for(State i = 0; i < *statesCnt; i++)
			buffer.emplace_back(i, LabelType::epsilon(), i);
		return *this;
	}
	std::span<const Transition<LabelType>> operator()(State st) const
	{
		if(!isSorted)
			throw std::runtime_error("cannot create subrange: transitions not sorted");
		return {buffer.begin() + startInd[st], buffer.begin() + startInd[st + 1]};
	}
	void clear()
	{
		startInd.clear();
		buffer.clear();
		isSorted = false;
	}
	friend std::ostream& operator<<(std::ostream& os, const TransitionList& list)
	{
		os << list.buffer.size() << '\n';
		for(auto&& tr : list.buffer)
			os << tr << '\n';
		os << list.isSorted << '\n';
		if(list.isSorted)
		{
			os << list.startInd.size() << '\n';
			for(std::size_t ind : list.startInd)
				os << ind << ' ';
			os << '\n';
		}
		return os;
	}
	friend std::istream& operator>>(std::istream& is, TransitionList& list)
	{
		list.clear();
		std::size_t transitionsCnt;
		is >> transitionsCnt;
		list.buffer.reserve(transitionsCnt);
		while(transitionsCnt--)
		{
			Transition<LabelType> tr;
			is >> tr;
			list.buffer.push_back(std::move(tr));
		}
		is >> list.isSorted;
		if(list.isSorted)
		{
			std::size_t startIndSize;
			is >> startIndSize;
			list.startInd.reserve(startIndSize);
			while(startIndSize--)
			{
				std::size_t ind;
				is >> ind;
				list.startInd.push_back(ind);
			}
		}
		return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
private:
	friend void sortByLabel<>(TransitionList<LabelType>& list);

	template<std::invocable<Transition<LabelType>> Proj>
	void sort(std::size_t maxValue, Proj proj)
	{
		isSorted = false;
		count(maxValue, proj);
		for(std::size_t i = 1; i < startInd.size(); i++)
			startInd[i] += startInd[i - 1];
		std::vector<Transition<LabelType>> sorted(buffer.size());
		for(auto& tr : buffer | std::views::reverse)
			sorted[--startInd[std::invoke(proj, tr)]] = std::move(tr);
		buffer = std::move(sorted);
	}
	void count(std::size_t maxValue, std::invocable<Transition<LabelType>> auto proj)
	{
		startInd.resize(maxValue + 2);
		std::ranges::fill(startInd, 0);
		for(const auto& tr : buffer)
			startInd[std::invoke(proj, tr)]++;
	}
};

template<class LabelType>
void sortByLabel(TransitionList<LabelType>& list)
{
	list.isSorted = false;
	std::ranges::sort(list.buffer, {}, &Transition<LabelType>::Label);
}

template<>
inline void sortByLabel(TransitionList<SymbolOrEpsilon>& list)
{
	list.sort(std::numeric_limits<std::make_unsigned_t<Symbol>>::max(), &Transition<SymbolOrEpsilon>::Label);
}

template<>
inline void sortByLabel(TransitionList<SymbolPair>& list)
{
	auto trLabelSecond = [](const Transition<SymbolPair>& tr) { return tr.Label().coords[1]; };
	list.sort(std::numeric_limits<std::make_unsigned_t<Symbol>>::max(), trLabelSecond);
	auto trLabelFirst = [](const Transition<SymbolPair>& tr) { return tr.Label().coords[0]; };
	list.sort(std::numeric_limits<std::make_unsigned_t<Symbol>>::max(), trLabelFirst);
}

template<>
inline void sortByLabel(TransitionList<Symbol_Word>& list)
{
	auto trLabelFirst = [](const Transition<Symbol_Word>& tr) { return tr.Label().first; };
	list.sort(std::numeric_limits<std::make_unsigned_t<Symbol>>::max(), trLabelFirst);
}

#endif
