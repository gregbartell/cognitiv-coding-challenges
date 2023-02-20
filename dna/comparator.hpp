#pragma once

#include "person.hpp"

#include <vector>
#include <utility>

namespace dna
{
	struct Difference
	{
		// [start, end) index of interesting segment
		using subsection = std::pair<size_t, size_t>;

		Difference(size_t chromosome_idx, size_t person_a_start, size_t person_a_end, size_t person_b_start, size_t person_b_end) : chromosome_idx{chromosome_idx}, person_a{person_a_start, person_a_end}, person_b{person_b_start, person_b_end}
		{}

		size_t chromosome_idx;
		subsection person_a;
		subsection person_b;
	};

	template <typename STREAM>
	STREAM& operator<<(STREAM& os, const Difference& d)
	{
		os << "Chromosome " << d.chromosome_idx;
		os << " | first sample: [" << d.person_a.first << ", " << d.person_a.second << "]";
		os << " second sample: [" << d.person_b.first << ", " << d.person_b.second << "]";
		return os;
	}

	class Comparator
	{
	private:
		// # of chromosomes in a valid sample
		static constexpr size_t NUM_CHROMOSOMES = 23;

		static constexpr size_t SEX_CHROMOSOME_IDX = 22;
		// Approx. lengths (in base pairs) of chromosome 23 for X/Y, used to determine genetic sex
		static constexpr size_t X_CHROMOSOME_LEN = 156'000'000;
		static constexpr size_t Y_CHROMOSOME_LEN =  57'000'000;

		// Fixed sequence of repeating bases in telomeres
		static constexpr std::array<base, 6> TELOMERE_SEQ = {T, T, A, G, G, G};

	public:
		enum class SexChromosome
		{
			X,
			Y,
			MAX, // Sentinel out-of-bounds value
		};

		// Result only valid when run on chromosome 23
		template <HelixStream H>
		static SexChromosome getSex(const H& helix)
		{
			if (auto helix_len = helix.size() * packed_size::value;
					(helix_len > (4*X_CHROMOSOME_LEN/5)) &&
					(helix_len < (5*X_CHROMOSOME_LEN/4)))
			{
				return SexChromosome::X;
			}
			else if (helix_len > (4*Y_CHROMOSOME_LEN/5) &&
					helix_len < (5*Y_CHROMOSOME_LEN/4))
			{
				return SexChromosome::Y;
			}

			// Bad length, unlikely to be a valid chromosome
			return SexChromosome::MAX;
		}

		// Returns [start, end) of the interesting data in a HelixStream
		// i.e. the data between telomeres
		// Returned values are indices of *bases*, NOT bytes
		template <HelixStream H>
		static std::pair<size_t, size_t> getDataRange(H& helix)
		{
			size_t telomere_idx = 0;

			size_t data_start = 0;
			size_t data_end = helix.size() * packed_size::value;

			// We need to find at least one complete telomere in order to classify it,
			// so if the helix is too short, we return early
			if (data_end < TELOMERE_SEQ.size())
			{
				return {data_start, data_end};
			}

			// How this function works:
			// 1. Identify partial telomere if present at beginning of data
			// 2. Adavance through data until it stops matching telomere pattern
			// 3. Repeat steps 1 and 2 in reverse for the end of the data

			helix.seek(0);
			auto buffer = helix.read();
			// XXX TODO
			// This function as written depends on the buffer being large enough to
			// hold the entire chromosome. This would clearly not be the case in real-world usage.
			// Adapting this code to handle streaming reads going forward is relatively simple,
			// but seeking backwards would be tougher. Perhaps a different approach altogether
			// would be better, but I've spent enough time on this already that I've decided to
			// submit what I've got so far.

			// For each possible starting position in the telomere sequence
			for (size_t t_idx = 0; t_idx < TELOMERE_SEQ.size(); t_idx++)
			{
				bool match = true;
				// Check each letter in the buffer against the corresponding position in the telomere sequence
				for (size_t b_idx = 0; b_idx < TELOMERE_SEQ.size() && b_idx < buffer.size(); b_idx++)
				{
					if (buffer[b_idx] != TELOMERE_SEQ[(t_idx + b_idx) % TELOMERE_SEQ.size()])
					{
						match = false;
						break;
					}
				}

				if (match)
				{
					telomere_idx = t_idx;
					data_start = TELOMERE_SEQ.size();
					break;
				}
			}

			// Keep iterating through the data until it stops matching telomeres
			while (data_start < data_end && buffer[data_start] == TELOMERE_SEQ[telomere_idx])
			{
				data_start++;
				telomere_idx = (telomere_idx + 1) % TELOMERE_SEQ.size();
			}

			// If there isn't enough room for a complete telomere at the end
			if (data_end < data_start + TELOMERE_SEQ.size())
			{
				return {data_start, data_end};
			}

			// Find the end of telomeres the same way as at the front, but backwards
			for (size_t t_idx = 0; t_idx < TELOMERE_SEQ.size(); t_idx++)
			{
				bool match = true;
				for (size_t b_idx = 0; b_idx < TELOMERE_SEQ.size(); b_idx++)
				{
					if (buffer[data_end - 1 - b_idx] != TELOMERE_SEQ[(TELOMERE_SEQ.size() + t_idx - b_idx) % TELOMERE_SEQ.size()])
					{
						match = false;
						break;
					}
				}

				if (match)
				{
					telomere_idx = t_idx;
					data_end -= TELOMERE_SEQ.size();
					break;
				}
			}

			// Keep iterating through the data until it stops matching telomeres
			while (data_end > data_start && buffer[data_end - 1] == TELOMERE_SEQ[telomere_idx])
			{
				data_end--;
				telomere_idx = (TELOMERE_SEQ.size() + telomere_idx - 1) % TELOMERE_SEQ.size();
			}

			return {data_start, data_end};
		}

		Comparator() = delete; // Static methods only, no instances should be constructed

		template <Person P>
		static std::vector<Difference> compare(const P& a, const P& b)
		{
			if (a.chromosomes() != NUM_CHROMOSOMES || b.chromosomes() != NUM_CHROMOSOMES)
			{
				throw std::invalid_argument("chromosome data does not match expected size");
			}

			std::vector<Difference> ret{};

			for (size_t chromosome_idx = 0; chromosome_idx < NUM_CHROMOSOMES; chromosome_idx++)
			{
				auto helix_a = a.chromosome(chromosome_idx);
				auto helix_b = b.chromosome(chromosome_idx);

				if (chromosome_idx == SEX_CHROMOSOME_IDX)
				{
					auto a_sex = getSex(helix_a);
					auto b_sex = getSex(helix_b);

					// We don't want to compare sex chromosomes if sexes are different
					if ((a_sex != b_sex) || (a_sex == SexChromosome::MAX))
					{
						continue;
					}
				}

				auto [a_start, a_end] = getDataRange(helix_a);
				auto [b_start, b_end] = getDataRange(helix_b);

				// XXX TODO: the general idea from here is to seek each stream to the start of the
				// interesting (non-telomere data) and compare chunks. With 99.9% of the genome
				// being the same between people, hopefully many chunks will be base-for-base
				// identical and not require any additional processing. Of the chunks that do,
				// I would look at something like this:
				// https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
				// to try to find large sequences that don't align well and report those
				// as meaningful differences between people.
			}

			return ret;
		}
	};
}
