#include "catch.hpp"
#include "fake_person.hpp"

#include "comparator.hpp"

struct FakeHelixStream
{
	FakeHelixStream(size_t sz) : m_size{sz} {}

	void seek(size_t) const {}
	void read() const {}
	auto size() const { return m_size; }

	size_t m_size = 0;
};

TEST_CASE("getSex classifies sexes")
{
	SECTION("Invalid length")
	{
		FakeHelixStream helix{0};
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::MAX);
		helix = 100'000'000;
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::MAX);
	}

	SECTION("X Chromosome")
	{
		FakeHelixStream helix{156'000'000 / dna::packed_size::value};
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::X);
		helix = (150'000'000 / dna::packed_size::value);
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::X);
		helix = (160'000'000 / dna::packed_size::value);
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::X);
	}

	SECTION("Y Chromosome")
	{
		FakeHelixStream helix{57'000'000 / dna::packed_size::value};
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::Y);
		helix = (50'000'000 / dna::packed_size::value);
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::Y);
		helix = (60'000'000 / dna::packed_size::value);
		CHECK(dna::Comparator::getSex(helix) == dna::Comparator::SexChromosome::Y);
	}
}

TEST_CASE("getDataRange identifies telomeres")
{
	using dna::A;
	using dna::C;
	using dna::G;
	using dna::T;

	// Note: most tests use repeated Cs as non-telomere data and other bases are part of telomeres

	std::vector<std::byte> data{};

	SECTION("Trivial (empty helix)")
	{
		fake_stream helix({}, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 0);
		CHECK(end == 0);
	}

	SECTION("No telomeres")
	{
		data = {dna::pack(C, C, C, C), dna::pack(C, C, C, C)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 0);
		CHECK(end == data.size() * dna::packed_size::value);
	}

	SECTION("Complete telomere at start")
	{
		data = {dna::pack(T, T, A, G), dna::pack(G, G, C, C)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 6);
		CHECK(end == data.size() * dna::packed_size::value);
	}

	SECTION("Multiple complete telomeres at start")
	{
		data = {dna::pack(T, T, A, G), dna::pack(G, G, T, T), dna::pack(A, G, G, G), dna::pack(C, C, C, C)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 12);
		CHECK(end == data.size() * dna::packed_size::value);
	}

	SECTION("Partial telomere at start")
	{
		data = {
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(C, C, C, C)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 20);
		CHECK(end == data.size() * dna::packed_size::value);
	}

	SECTION("Complete telomere at end")
	{
		data = {
			dna::pack(C, C, C, C),
			dna::pack(C, C, T, T),
			dna::pack(A, G, G, G)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 0);
		CHECK(end == 6);
	}

	SECTION("Multiple complete telomeres at end")
	{
		data = {dna::pack(C, C, C, C), dna::pack(T, T, A, G), dna::pack(G, G, T, T), dna::pack(A, G, G, G)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 0);
		CHECK(end == 4);
	}

	SECTION("Partial telomere at end")
	{
		data = {
			dna::pack(C, C, C, C),
			dna::pack(C, C, C, C),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 0);
		CHECK(end == 8);
	}

	SECTION("Partial telomeres at start and end")
	{
		data = {
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(C, C, C, C),
			dna::pack(C, C, C, C),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 20);
		CHECK(end == 28);
	}

	// This test is the same as the above, but the non-telomere data uses T/A/G to
	// make sure they don't get misinterpreted as telomeres
	SECTION("Telomere-like data between telomeres")
	{
		data = {
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(G, G, G, G),
			dna::pack(T, T, T, T),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T)};
		fake_stream helix(data, 128);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 20);
		CHECK(end == 28);
	}

	// XXX TODO: this test may or may not pass, but this functionality is not properly
	// implemented yet and this code reads out-of-bounds memory

	// This test is designed to be sure that we can still process data even when using a
	// buffer chunk size too small to hold an entire telomere
	SECTION("Partial telomeres at start and end; small helix buffer chunk size")
	{
		data = {
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T),
			dna::pack(A, G, G, G),
			dna::pack(C, C, C, C),
			dna::pack(C, C, C, C),
			dna::pack(C, C, T, T),
			dna::pack(A, G, G, G),
			dna::pack(T, T, A, G),
			dna::pack(G, G, T, T)};
		fake_stream helix(data, 1);

		auto [start, end] = dna::Comparator::getDataRange(helix);
		CHECK(start == 20);
		CHECK(end == 30);
	}
}
