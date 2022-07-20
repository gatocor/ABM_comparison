/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void setup_microenvironment( void )
{
	// set domain parameters
	
	initialize_microenvironment(); 	
	
	// these will ***overwrite*** values specified in the 
	// microenvironment_setup part of the XML,
	// based on what's in the user_parameters section 
	
	microenvironment.name = "synthetic tissue"; 
	
	// int cargo_index = microenvironment.find_density_index( "cargo signal" ); 
	// int director_index = microenvironment.find_density_index( "director signal" ); 
	
	// microenvironment.diffusion_coefficients[cargo_index] = 
	// 	parameters.doubles("cargo_signal_D");  
	// microenvironment.decay_rates[cargo_index] = 
	// 	parameters.doubles("cargo_signal_decay");  
	
	// microenvironment.diffusion_coefficients[director_index] = 
	// 	parameters.doubles("director_signal_D");  
	// microenvironment.decay_rates[director_index] = 
	// 	parameters.doubles("director_signal_decay"); 
	
	// display the microenvironment again 
	
	microenvironment.display_information( std::cout ); 
	
	return; 
}

void create_cell_types( void )
{
	SeedRandom( parameters.ints("random_seed") ); 
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// turn off proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0; 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 
	
	// int cargo_index = microenvironment.find_density_index( "cargo signal" ); // 1 
	// int director_index = microenvironment.find_density_index( "director signal" ); // 0 
	
	// set uptake and secretion to zero 
	// cell_defaults.phenotype.secretion.secretion_rates[director_index] = 0; 
	// cell_defaults.phenotype.secretion.uptake_rates[director_index] = 0; 
	// cell_defaults.phenotype.secretion.saturation_densities[director_index] = 1; 
	
	// cell_defaults.phenotype.secretion.secretion_rates[cargo_index] = 0; 
	// cell_defaults.phenotype.secretion.uptake_rates[cargo_index] = 0; 
	// cell_defaults.phenotype.secretion.saturation_densities[cargo_index] = 1; 

	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = NULL; 
		
	//
	// Define "seed" cells 
	
	cell = cell_defaults; 
	cell.type = cell_ID; 
	cell.name = "cell"; 
		
	// seed cell rule 
	
	cell.functions.update_phenotype = NULL; 
		
	cell.phenotype.motility.is_motile = true; 
	cell.phenotype.motility.persistence_time = 
		parameters.doubles("motility_persistence_time"); // 5.0; 
	cell.phenotype.motility.migration_speed = 
		parameters.doubles("migration_speed"); // 5; 
	cell.phenotype.motility.migration_bias = 
		parameters.doubles("migration_bias"); // 5.0; 

	cell.functions.update_migration_bias = NULL;

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

std::vector<std::string> robot_coloring_function( Cell* pCell )
{
	std::string color = "black"; 
	std::vector< std::string > output( 4 , color ); 
	
	// black cells if necrotic 
	if( pCell->phenotype.death.dead == true )
	{ return output; }

	output[3] = "none"; // no nuclear outline color 
	
	static std::string cell_color = parameters.strings( "cell_color" ); 

	if( pCell->type == cell_ID )
	{ color = cell_color; }
	
	output[0] = color; 
	output[2] = color; 
	
	return output; 
}

void setup_tissue( void )
{
	int number_of_cells = parameters.ints("number_of_cells"); // 15;  

	std::cout << "Placing cells ... " << std::endl; 
	
	// randomly place seed cells 
	
	std::vector<double> position(3,0.0); 
	
	double x_range = default_microenvironment_options.X_range[1] - default_microenvironment_options.X_range[0]; 
	double y_range = default_microenvironment_options.Y_range[1] - default_microenvironment_options.Y_range[0]; 

	double relative_margin = 0.2;  
	double relative_outer_margin = 0.02; 
	
	for( int i=0; i < number_of_cells ; i++ )
	{
		// pick a random location 
		position[0] = default_microenvironment_options.X_range[0] + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = default_microenvironment_options.Y_range[0] + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 
		
		// place the cell
		Cell* pC;
		pC = create_cell( cell ); 
		pC->assign_position( position );
		pC->is_movable = true; 
	}	

	std::cout << "done!" << std::endl; 
	// make a plot 
	
	PhysiCell_SVG_options.length_bar = 200; 
	SVG_plot( "initial.svg" , microenvironment, 0.0 , 0.0 , robot_coloring_function );	
	
	return; 
}

// void cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
// {
// 	// if attached, biased motility towards director chemoattractant 
// 	// otherwise, biased motility towards cargo chemoattractant 
	
// 	static double attached_migration_bias = 
// 		parameters.doubles("attached_migration_bias"); 
// 	static double unattached_migration_bias = 
// 		parameters.doubles("unattached_migration_bias"); 
		
// 	static int cargo_index = microenvironment.find_density_index( "cargo signal" ); // 1 
// 	static int director_index = microenvironment.find_density_index( "director signal" ); // 0 
	
// 	if( pCell->state.neighbors.size() > 0 )
// 	{
// 		phenotype.motility.migration_bias = attached_migration_bias; 

// 		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(director_index);	
// 		normalize( &( phenotype.motility.migration_bias_direction ) );			
// 	}
// 	else
// 	{
// 		phenotype.motility.migration_bias = unattached_migration_bias; 
		
// 		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(cargo_index);	
// 		normalize( &( phenotype.motility.migration_bias_direction ) );			
// 	}
	
// 	return; 
// }
