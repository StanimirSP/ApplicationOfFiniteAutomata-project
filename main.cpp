#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include "regularExpression.h"
#include "ThompsonsConstruction.h"
#include "transducer.h"
#include "replacement.h"
#include "bimachine.h"

std::string readFromFile(const std::filesystem::path& path, char delim = '\n')
{
	std::ifstream ifs(path);
	if(!ifs)
		throw std::runtime_error("could not open \"" + path.string() + "\" for reading");
	std::string str;
	std::getline(ifs, str, delim);
	return str;
}

int main(int argc, char** argv) try
{
	using namespace std::string_literals;
	if(argc != 4)
		throw std::runtime_error("Usage: "s + (argv[0] ? argv[0] : "<program name>") + " <alphabet_file> <pattern_file> <text_file>");
	std::string alphabet = readFromFile(argv[1], Constants::Epsilon);
	Transducer<false, Symbol_Word> rtime;
	{
		std::cerr << "Constructing regular expression...\n";
		RegularExpression<WordPair> r(readFromFile(argv[2]));
		std::cerr << "Constructing letter transducer...\n";
		LetterTransducer lt = Transducer(regexToMFSA(std::move(r), alphabet)).expand();
		std::cerr << "Constructing transducer for LML...\n";
		rtime = ReplaceLML(lt, alphabet).realTime();
	}
	std::cerr << "Constructing bimachine...\n";
	Bimachine bm(std::move(rtime), alphabet);
	//std::cerr << bm << '\n';
	std::string text = readFromFile(argv[3], '\0');
	std::cerr << "Replacing text...\n";
	std::cout << bm(text) << std::endl;

	/*{
		//RegularExpression<WordPair> r("([aa,x]|[a,yy])([aa,x]|[a,yy])*");
		//RegularExpression<WordPair> r("[a,_][_,x]*"); // infinitely ambiguous
		//RegularExpression<WordPair> r("[a,x]|[a,y]|[a,z]|[x,x][a,t][a,t]*");
		//RegularExpression<WordPair> r("[ ,_][ ,_]*[_, ]|[a,A]|[b,B]|[c,C]|[d,D]|[e,E]|[f,F]|[new york university,01]|[new york,02]|[york,03]|[york university,044]");
		//RegularExpression<WordPair> r("[`,_][`,_]*[_,`]|[a,A]|[b,B]|[c,C]|[d,D]|[e,E]|[f,F]");
		//RegularExpression<WordPair> r("[`,_][`,_]*[_,`]|[a,A]|[b,B]|[d,D]");
		//RegularExpression<WordPair> r("[x,y][x,y]*[y,x]|[a,A]|[b,B]|[c,C]|[d,D]|[e,E]|[f,F]");
		std::string s;
		std::getline(std::cin, s);
		//RegularExpression<WordPair> r(s);
		std::string alphabet("abcfdefghijklmn`o pqrstuvw!yz.x0123456789");
		//std::string alphabet("abd`");
		//std::cerr << alphabet << '\n';
		LetterTransducer T3(Transducer(regexToMFSA(std::move(r), alphabet)).expand());
		//std::cout << T3.pseudoMinimize() << '\n';
		//using namespace std::string_view_literals;
		//auto R = Replace(T3, "abyx"sv);
		//R.trim().transitions.sort();
		//std::cout << R << '\n';
		auto lml = ReplaceLML(T3, alphabet);
		lml.pseudoMinimize();
		//lml.transitions.sort();
		//std::cout << lml << '\n';
		//std::cout << lml.pseudoMinimize() << '\n';
		auto rt = lml.realTime();
		//rt.transitions.sort();
		//std::cout << rt << '\n';
		Bimachine bm(std::move(rt), alphabet);

		std::cout << bm(s) << std::endl;
		//std::cout << bm("xaaaaxy") << std::endl;
		//std::cout << bm("sd") << std::endl;
	}*/
	/*{
		RegularExpression<SymbolOrEpsilon> r("a*bb*c"), r2("c*|aabbc*");
		std::cout << r.TokenizedRegex() << '\n';
		std::cout << r.TokenizedReversePolishNotation() << '\n';
		for(auto&& x : r.BaseTokens())
			std::cout << x << '\n';
		ClassicalFSA ca = regexToMFSA(std::move(r), "abc"), ca2 = regexToMFSA(std::move(r2), "abc");
		std::cout << ca << '\n';
		//std::cout << ca.convertToDFSA() << '\n';
		//std::cout << ca2.convertToDFSA() << '\n';
		std::cout << ca.intersect(ca2) << '\n';
	}*/
	/*{
		RegularExpression<SymbolOrEpsilon> r("(a|b)ba*b");

		ClassicalFSA T = regexToMFSA(std::move(r), "ab");
		//T.convertToDFSA();
		std::cerr << T << '\n';
		std::cerr << T.pseudoMinimize() << '\n';
	}*/

	/*{
		RegularExpression<SymbolOrEpsilon> r("(a|b)*a(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b|_)(a|b|_)(a|b|_)(a|b|_)(a|b|_)(a|b|_)");
		RegularExpression<SymbolOrEpsilon> r2("(a|b)*a(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)(a|b)");

		ClassicalFSA T = regexToMFSA(std::move(r), "ab");
		ClassicalFSA T2 = regexToMFSA(std::move(r2), "ab");
		ClassicalFSA result = T.intersect(T2.complement());
		std::cerr << "Minimizing...\n";
		result.pseudoMinimize().print();
	}*/
}
catch(const std::exception& e)
{
	std::cerr << e.what() << '\n';
	return 1;
}
