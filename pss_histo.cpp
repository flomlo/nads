#include <fstream>
#include <iostream>
#include <cstdint>

#define USE_CELLS_WITHOUT_DIMENSION

#include "include/argparser.h"
#include "include/parameters.h"
#include "include/usage/flagser-count.h"
#include "include/complex/directed_flag_complex_in_memory_computer.h" 
//#include "include/complex/directed_flag_complex_computer.h" 


typedef std::vector<std::pair<int, vertex_index_t>> coboundaries_t;
typedef struct{} NoneType;

struct pss_computer_t {
  const filtered_directed_graph_t& graph;
  uint64_t pss = 0;
  coboundaries_t coboundaries{};  //have on coboundary vector per thread, reuse per simplex to save allocs

  std::vector<size_t> pss_histo; //{std::vector<size_t>(10,0)}; //histogram of postsynaptic d+1 simplices per d-simplex 

  void done() {}

  void operator()(vertex_index_t* first_vertex, int size, NoneType) {
    std::vector<size_t> vertex_offsets;
    for (int j = 0; j < size; j++) vertex_offsets.push_back(first_vertex[j] >> 6);

    coboundaries.clear();
    
    // find all coboundaries the the simplex given by 
    for (int i = 0; i <= size; i++) {
      // Check intersections in chunks of 64
      for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
        size_t bits = -1; // All bits set

        for (int j = 0; bits > 0 && j < size; j++) {
          // Remove the vertices already making up the cell
          if (vertex_offsets[j] == offset) bits &= ~(ONE_ << (first_vertex[j] - (vertex_offsets[j] << 6)));

          // Intersect with the outgoing/incoming edges of the current vertex
          bits &= j < i ? graph.get_outgoing_chunk(first_vertex[j], offset)
                        : graph.get_incoming_chunk(first_vertex[j], offset);
        }

        size_t vertex_offset = offset << 6;
        while (bits > 0) {
          // Get the least significant non-zero bit
          auto b = __builtin_ctzl(bits);

          // Unset this bit
          bits &= ~(ONE_ << b);

          coboundaries.push_back({i,vertex_offset+b});
        }
      }
    }
    // now calculate the number of fully postsynaptic simplices
    uint64_t new_pss = 0;
    for (auto cb : coboundaries) {
      if (cb.first == size) new_pss += 1;
    }
    pss += new_pss;
    // add value to histogram
    if(new_pss >= pss_histo.size()) pss_histo.resize(new_pss+1);
    pss_histo[new_pss] += 1;
  }
};


std::pair<std::vector<uint64_t>, std::vector<std::vector<size_t>>> count_pss(const filtered_directed_graph_t& graph, const flagser_parameters& params) {
  int min_dimension = 1;
  if (params.min_dimension >= 1) min_dimension = params.min_dimension;
  int max_dimension = 100; //hacky
  if (params.max_dimension >= 1) max_dimension = params.max_dimension;
  assert(min_dimension <= max_dimension);

  // generate flag complex
  directed_flag_complex_in_memory_t<NoneType> complex(graph, params.nb_threads, max_dimension);
  
  std::cout << "generated flag complex" << std::endl;

	std::vector<uint64_t> pss;
  std::vector<std::vector<size_t>> pss_histo;
  for (int dim = min_dimension-1; dim <= max_dimension-1; dim++) {
    //initialize threads
		std::vector<pss_computer_t> pss_counter(params.nb_threads, pss_computer_t{graph});
    //actually compute
    complex.for_each_cell(pss_counter, dim); 
    
    //sum results
    uint64_t pss_in_dim = 0;
    for (const auto& thread : pss_counter) pss_in_dim += thread.pss;

    //sum histograms
    size_t hist_max_size = 0;
    for (const auto& thread : pss_counter) hist_max_size = std::max(thread.pss_histo.size(), hist_max_size);
    std::vector<size_t> pss_histo_in_dim(hist_max_size, 0);
    for (const auto& thread : pss_counter) {
      for (size_t i = 0; i < hist_max_size; i++) {
        if (i < thread.pss_histo.size()) pss_histo_in_dim[i] += thread.pss_histo[i];
      }
    }
    // pss_in_dim == 0 => pss in higher dimensions must also be 0, thus we may abort.
    if (pss_in_dim == 0) break;
    
    //debug
    /*
    std::cout << '[';
    for (size_t i = 0; i < hist_max_size; i++) {
      std::cout << pss_histo_in_dim[i] << ',' << ' ';
    }
    std::cout << ']' ;
    */

    std::cout << dim+1 << ": " << pss_in_dim << std::endl;

    pss.push_back(pss_in_dim);
    pss_histo.push_back(pss_histo_in_dim);
  }
	return std::make_pair(pss, pss_histo);
}

int main(int argc, char** argv) {
	try {
		auto arguments = parse_arguments(argc, argv);

		auto positional_arguments = get_positional_arguments(arguments);
		auto named_arguments = get_named_arguments(arguments);
		auto params = flagser_parameters(named_arguments);
		named_arguments_t::const_iterator it;

		const char* input_filename = positional_arguments[0];

		filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, params);

		auto cell_count = count_pss(graph, params);
	} catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}
