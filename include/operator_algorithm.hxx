//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <rsa_ghost_area.hxx>
#include <rsa_domain.hxx>
#include <rsa_random.hxx>
#include <list_of_voxels.hxx>
#include <operator_paraview.hxx>
#include <operator_check_cell.hxx>
#include <operator_ghost.hxx>
#include <radius_generator.hxx>

namespace algorithm {
using namespace sac_de_billes;
using namespace sac_de_billes;

template<int DIM>
class rsa_algo {
	//! @brief Class for execution of the algorithm.
public:
	//! @brief class for proceeding with the rsa_algo
	//! @param a_domain : contains all the cells and spheres
	//! @param a_radius_generator : a functor for producing new radii
	//! @param a_nbshots_singledraw : nb of local draws at each mpi proc, before exchanging with other procs
	//! @param a_n_draw : for naive method, nb of draws
	rsa_algo(rsa_domain<DIM>& a_domain, RadiusGenerator<DIM>& a_radius_generator,
		int a_nbshots_singledraw, int a_n_draw = 10) :
		m_nbshots_singledraw{ std::move(a_nbshots_singledraw) },
		m_n_draw{ a_n_draw },
		m_radius_generator{ a_radius_generator },
		m_domain{ a_domain },
		m_ghost_data{},
		m_r_max{ a_radius_generator.get_max_radius() } {
		this->get_grid() = rsa_grid<DIM>(a_radius_generator.get_max_radius(), a_domain.get_ghost_layer(),
			a_domain.get_inf(), a_domain.get_sup());
		if (a_radius_generator.get_max_radius() > a_domain.get_m_rad()) {
			throw runtime_error("Impossible to have the maximal radius of radius generator larger than domain implicit radius");
		}
	}

private:
	//! @brief desired number of shots for a single draw (idea : 1 draw = 1 MPI communication)
	int m_nbshots_singledraw;
	//! @brief nb of draw. Useful only for the naive method.
	int m_n_draw;
	//! @brief radius generator
	RadiusGenerator<DIM> m_radius_generator;
	//! @brief domain in which the spheres are drawn
	rsa_domain<DIM>& m_domain;
	//! @brief storage for ghost data
	rsa_data_storage<DIM> m_ghost_data;
	//! @brief r_max
	double m_r_max;

public:
	//! @brief Execute the RSA algorithm.
	//! @tparam method 1 is standard. 0 is for the naive method.
	//! @param seed :for the pseudo-random generator.
	template<int method>
	void proceed(size_t seed);
	//! @brief: Reset the radius generator, but retains the already placed spheres.\n
	//! Proceed can then be appealed once more.
	void reset_radius_generator(RadiusGenerator<DIM>& a_radius_generator);

	//! @brief Gets the grid of cells.
	//! @return The rsa_grid<DIM> representing the grid of cells.
	const rsa_grid<DIM>& get_grid() const { return m_domain.get_grid(); }
	//! @brief Gets the grid of cells.
	//! @return The rsa_grid<DIM> representing the grid of cells.
	rsa_grid<DIM>& get_grid() { return m_domain.get_grid(); }

	//! given an already packed configuration, checks how much time the unavoidable voxellation process only requires
	void check_vox_time();

private:
	//! @brief naive method without an underlying grid for performance comparison
	//! @param random_generator 
	void proceed_naive(std::mt19937& random_generator);
	//! @brief best method from Ebeida et al
	//! @param random_generator 
	void proceed_voxel(std::mt19937& random_generator);

	//! @brief  Comput the (average) number of spheres that should be drawn, based on the number of voxels
	//! @param number_of_voxels : how many free voxels are there?
	//! @return : average number of sphere that should be drawn.
	double compute_intensity_poisson(int number_of_voxels);
	//! compute the desired number of shots
	//! renormalizes it, in such a way that, the (expectance) of the maximal number of shots is < m_desired_nb_draw
	//! @see  m_desired_nb_draw
	//! @param intensity : desired intensity
	//! @param maximum_intensity : maximal intensity accross all the MPI processes
	int compute_nb_shots_voxel(std::mt19937& random_generator, double intensity, double maximum_intensity);
	//! @brief decide to remove covered voxels or subdivide them, depending on miss_rate
	//! @param miss_rate : observed miss rate in shots
	//! @param uncovered_voxels : uncovered voxels for drawing random spheres
	void update_covered_voxels(voxel_list::list_of_voxels<DIM>& uncovered_voxels,
		double miss_rate);
	//! @see : compute_intensity_poisson
	static constexpr double  m_intensity_factor = 0.01; // magical constant
	//static constexpr double  m_intensity_factor = 0.05; // magical constant
	//! @brief : to do
	//! @return : nb of newly drawn spheres
	//! @param center_generator : generates the center of spheres
	//! @param priority_generator : generate the priority of spheres, to know which one come first
	//! @param nb_shots : desired nb of centers to be drawn, for each MPI process
	//! @param may_outreach_nb_spheres : true if the number of shots is higher than the desired nb of spheres to be placed
	template<class CenterGenerator, class PriorityGenerator>
	int64_t single_draw(CenterGenerator& center_generator, PriorityGenerator& priority_generator, int nb_shots,
		bool may_outreach_nb_spheres = true);

	//! @brief update the radius generator
	void update_radius_generator(int64_t nb_placed_spheres);
	//! @return return the desired miss rate of the voxel strategy. (Magical constants.)
	double desired_miss_rate() const;
};

//! @brief Adds a sample to the main grid.
//! @param a_grid : The main grid for data storage.
//! @param a_recv : The buffers used to store received data from neighboring processes.
//! @param a_send : The buffers used store data to be sent to neighboring processes.
//! @param a_ghost_areas : The ghost areas used to store ghost data.
//! @param a_cells_with_conflicts : List of cells containing conflicts.
//! @return A boolean indicating if there is a conflicted spheres or not.
template<int DIM>
bool add_to_sample(
	rsa_grid<DIM>& a_grid,
	std::vector<buffer_for_spheres<DIM>>& a_recv, std::vector<buffer_for_spheres<DIM>>& a_send,
	std::vector<rsa_ghost_area<DIM>>& a_ghost_areas,
	std::vector<uint64_t>& a_cells_with_conflicts,
	uint64_t nb_spheres_max,
	int64_t& local_nb_spheres,
	bool may_outreach_nb_spheres);
//! @brief Generates spheres according to a center generator.
//! @return : Nb of newly added spheres
//! @param a_cells : The main grid of cells (validated spheres).
//! @param a_recv : The buffers storing received data from neighboring processes.
//! @param a_send : The buffers storing data to be sent to neighboring processes.
//! @param a_ghost_areas : The ghost areas used to store ghost.
//! @param a_rad : The radius of the spheres to be generated.
//! @param a_center_generator : The center generator used to create sphere centers.
//! @param a_size : The size of the sample.
//! @param a_ghost_data : The data storage used to store temporarly ghost spheres.
template<int DIM, typename CenterGenerator, class PriorityGenerator>
int64_t generate_spheres(rsa_grid<DIM>& a_grid,
	Buffers<DIM>& a_recv, Buffers<DIM>& a_send,
	GhostAreas<DIM>& a_ghost_areas,
	CenterGenerator& a_center_generator,
	RadiusGenerator<DIM>& a_radius_generator,
	PriorityGenerator& a_priority_generator,
	int a_size,
	rsa_data_storage<DIM>& a_ghost_data,
	bool may_outreach_nb_spheres);

//! @brief Generate candidate spheres
//! @tparam CenterGenerator
//! @tparam DIM
//! @param a_data
//! @param a_radius_generator
//! @param a_center_generator
//! @param a_size : nb of generated spheres
template<int DIM, typename CenterGenerator, class PriorityGenerator>
rsa_data_storage<DIM> generate_candidates(
	CenterGenerator& a_center_generator,
	const RadiusGenerator<DIM>& a_radius_generator,
	PriorityGenerator& a_priority_generator,
	int a_size);

//! @brief : This function generates a_size * a_n_draw spheres and try to add it in data_storage
template<int DIM, int method>
rsa_algo<DIM> uniform_generate(rsa_domain<DIM>& a_domain, double a_rad,
	int a_size, int a_n_draw = 1, size_t seed = 0);
template<int DIM, int method>
rsa_algo<DIM> uniform_generate(rsa_domain<DIM>& a_domain,
	sac_de_billes::RadiusGenerator<DIM>& radius_generator,
	int a_size, int a_n_draw, size_t seed = 0);

//! @brief : Only performs the unavoidable voxellation processes, with a final configuration.
template<int DIM>
void check_time_unavoidable_vox(rsa_domain<DIM>& a_domain, double min_radius, double exclusion_distance);

namespace auxi {
//! @brief generate candidate center positions for a given number of spheres
//! @tparam FuncLaw : type of double random generator
//! @tparam DIM : dimension of space
//! @param a_func_law : DIM-D array of random double generator for generating a position
//! @param a_size : nb of spheres position to generate
//! @return : a vector of sphere positions
template<int DIM, typename FuncLaw>
arr_vec_double<DIM> generate_sphere_positions(FuncLaw* a_func_law, int a_size);
//! @return a random law for generating a random position uniformly in the given domain
//! @tparam DIM : space dimension
//! @param a_domain : domain
template<int DIM>
std::vector<law::uniform<double>> create_random_law(const rsa_domain<DIM>& a_domain,
	std::mt19937& generator);
//! @brief only retain the first "a_maximal_nb_spheres" among the Undecided spheres.
//! @return the number of spheres that are selected.
template<int DIM>
int64_t retain_only_first_elements(rsa_grid<DIM>& a_grid,
	std::vector<uint64_t>& a_cells_with_conflicts,
	const uint64_t& a_maximal_nb_spheres);

//! update the sphere tags, putting  
//! - Conflict into Undecided
//!	- Undecided into Validated
//! - Not_Validated -> remove
//! - Skip -> remove
template<int DIM>
void update_sphere_tags(rsa_grid<DIM>& a_grid,
	std::vector<buffer_for_spheres<DIM>>& a_recv, std::vector<buffer_for_spheres<DIM>>& a_send,
	std::vector<rsa_ghost_area<DIM>>& a_ghost_areas,
	std::vector<uint64_t>& a_cells_with_conflicts,
	uint64_t nb_undecided);

//! @return : all the cell_id (relative to the grid) containing the center of a sphere (relative to spheres)
//! @param a_spheres
//! @param a_grid
template<int DIM>
vector<uint64_t> compute_conflict_cells(const rsa_data_storage<DIM>& a_spheres, const rsa_grid<DIM>& a_grid);
//! @return : all the cell_id (relative to the grid) containing a center of a sphere not validated
//! @param a_grid
template<int DIM>
void recompute_conflict_cells(const rsa_grid<DIM>& a_grid, std::vector<uint64_t>& a_cells_with_conflicts);
} // namespace auxi

//! @brief Generates a std::vector<INT_TYPE> of priorities.
//! @param a_size : The size of the std::vector<INT_TYPE> to be generated.
//! @param a_random_generator : random generator
//! @return The std::vector<INT_TYPE> containing the generated priorities.
template<class INT_TYPE>
std::vector <uint64_t> generate_priority(int a_size, std::mt19937& a_random_generator);
} // namespace algorithm

#include <operator_algorithm.ixx>
