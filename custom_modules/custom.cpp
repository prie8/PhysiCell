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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Global variables for substrate and cell type indices
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int glucose_substrate_index = -1;
int biofilm_cell_type_index = -1;
int probiotic_cell_type_index = -1;
int dnase_substrate_index = -1;
int il10_substrate_index = -1;
int prey_apoptosis_index = -1; // Index for the prey's apoptosis death model

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// create_cell_types function
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_cell_types( void )
{
  // set the random seed
  if (parameters.ints.find_index("random_seed") != -1) { SeedRandom(parameters.ints("random_seed")); }
  else { SeedRandom(); }

  initialize_default_cell_definition();
  cell_defaults.functions.volume_update_function = standard_volume_update_function;
  cell_defaults.functions.update_velocity = standard_update_cell_velocity;
  cell_defaults.functions.update_phenotype = NULL;
  cell_defaults.functions.custom_cell_rule = NULL;
  cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
  cell_defaults.functions.calculate_distance_to_membrane = NULL;

  initialize_cell_definitions_from_pugixml();

  cell_defaults.functions.update_phenotype = phenotype_function;
  cell_defaults.functions.custom_cell_rule = custom_function;

  Cell_Definition* pPreyDef = find_cell_definition( "S_aureus_biofilm" );
  Cell_Definition* pPredDef = find_cell_definition( "probiotic_E_coli" );

  bool prey_def_found = (pPreyDef != NULL);
  bool pred_def_found = (pPredDef != NULL);

  if(prey_def_found) {
    pPreyDef->functions.custom_cell_rule = prey_custom_function;
    pPreyDef->functions.update_phenotype = prey_phenotype_function;
    pPreyDef->functions.update_migration_bias = prey_motility_function;
    biofilm_cell_type_index = pPreyDef->type;
    prey_apoptosis_index = pPreyDef->phenotype.death.find_death_model_index("apoptosis");
    if (prey_apoptosis_index < 0) { std::cerr << "Warning: Apoptosis death model not found for S_aureus_biofilm!" << std::endl; }
  } else { std::cerr << "ERROR: Could not find Cell_Definition for S_aureus_biofilm!" << std::endl; }

  if(pred_def_found) {
    pPredDef->functions.custom_cell_rule = predator_custom_function;
    pPredDef->functions.update_phenotype = predator_phenotype_function;
    pPredDef->functions.update_migration_bias = predator_motility_function;
    probiotic_cell_type_index = pPredDef->type;
  } else { std::cerr << "ERROR: Could not find Cell_Definition for probiotic_E_coli!" << std::endl; }

  build_cell_definitions_maps();
  setup_signal_behavior_dictionaries();
  setup_cell_rules();
  display_cell_definitions( std::cout );

  return;
}
void create_cell_types1( void )
{
  // set the random seed
  if (parameters.ints.find_index("random_seed") != -1)
  {
    SeedRandom(parameters.ints("random_seed"));
  } else {
    SeedRandom(); // Seed with default if not specified
  }

  initialize_default_cell_definition();

  cell_defaults.functions.volume_update_function = standard_volume_update_function;
  cell_defaults.functions.update_velocity = standard_update_cell_velocity;
  cell_defaults.functions.update_phenotype = NULL;
  cell_defaults.functions.custom_cell_rule = NULL;
  cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
  cell_defaults.functions.calculate_distance_to_membrane = NULL;

  initialize_cell_definitions_from_pugixml();

  cell_defaults.functions.update_phenotype = phenotype_function;
  cell_defaults.functions.custom_cell_rule = custom_function;

  Cell_Definition* pPreyDef = find_cell_definition( "S_aureus_biofilm" );
  Cell_Definition* pPredDef = find_cell_definition( "probiotic_E_coli" );

  bool prey_def_found = (pPreyDef != NULL);
  bool pred_def_found = (pPredDef != NULL);

  // Assign functions and store indices for S_aureus_biofilm (Prey)
  if(prey_def_found)
  {
    pPreyDef->functions.custom_cell_rule = prey_custom_function;
    pPreyDef->functions.update_phenotype = prey_phenotype_function;
    pPreyDef->functions.update_migration_bias = prey_motility_function;
    biofilm_cell_type_index = pPreyDef->type; // Store index
                                              // *** Get the apoptosis model index for the prey ***
    prey_apoptosis_index = pPreyDef->phenotype.death.find_death_model_index("apoptosis");
    if (prey_apoptosis_index < 0) {
      std::cerr << "Warning: Apoptosis death model not found for S_aureus_biofilm!" << std::endl;
    }
  }
  else
  {
    std::cerr << "ERROR: Could not find Cell_Definition for S_aureus_biofilm!" << std::endl;
  }

  // Assign functions and store index for probiotic_E_coli (Predator)
  if(pred_def_found)
  {
    pPredDef->functions.custom_cell_rule = predator_custom_function;
    pPredDef->functions.update_phenotype = predator_phenotype_function;
    pPredDef->functions.update_migration_bias = predator_motility_function;
    probiotic_cell_type_index = pPredDef->type; // Store index
  }
  else
  {
    std::cerr << "ERROR: Could not find Cell_Definition for probiotic_E_coli!" << std::endl;
  }

  build_cell_definitions_maps();
  setup_signal_behavior_dictionaries();
  setup_cell_rules();
  display_cell_definitions( std::cout );

  return;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// setup_microenvironment function
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void setup_microenvironment( void )
{
  if( default_microenvironment_options.simulate_2D == true ) {
    std::cout << "Warning: overriding XML config options with 2D defaults!\n";
    default_microenvironment_options.X_range = {-500, 500};
    default_microenvironment_options.Y_range = {-500, 500};
    default_microenvironment_options.simulate_2D = true;
  }

  initialize_microenvironment();

  // *** Get substrate indices ***
  glucose_substrate_index = microenvironment.find_density_index("glucose");
  dnase_substrate_index = microenvironment.find_density_index("dnase");   // *** NEW ***
  il10_substrate_index = microenvironment.find_density_index("il10");    // *** NEW ***

  // Check if indices were found
  if (glucose_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'glucose'!" << std::endl; }
  if (dnase_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'dnase'!" << std::endl; }
  if (il10_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'il10'!" << std::endl; }

  return;
}
void setup_microenvironment1( void )
{
  if( default_microenvironment_options.simulate_2D == true )
  {
    std::cout << "Warning: overriding XML config options with 2D defaults!\n";
    default_microenvironment_options.X_range = {-500, 500};
    default_microenvironment_options.Y_range = {-500, 500};
    default_microenvironment_options.simulate_2D = true;
  }

  initialize_microenvironment();

  glucose_substrate_index = microenvironment.find_density_index("glucose");
  if (glucose_substrate_index < 0) {
    std::cerr << "ERROR: Could not find microenvironment density 'glucose'!" << std::endl;
  }

  return;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// setup_tissue function
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void setup_tissue( void )
{
  double Xmin = microenvironment.mesh.bounding_box[0];
  double Ymin = microenvironment.mesh.bounding_box[1];
  double Zmin = microenvironment.mesh.bounding_box[2];
  double Xmax = microenvironment.mesh.bounding_box[3];
  double Ymax = microenvironment.mesh.bounding_box[4];
  double Zmax = microenvironment.mesh.bounding_box[5];
  if( default_microenvironment_options.simulate_2D == true ) { Zmin = 0.0; Zmax = 0.0; }
  double Xrange = Xmax - Xmin;
  double Yrange = Ymax - Ymin;
  double Zrange = Zmax - Zmin;

  Cell* pC;

  // place prey (S_aureus_biofilm)
  Cell_Definition* pCD_Prey = find_cell_definition( "S_aureus_biofilm");
  if (!pCD_Prey) {
    std::cerr << "ERROR: Could not find Cell_Definition for S_aureus_biofilm during tissue setup!" << std::endl;
  } else {
    std::cout << "Placing cells of type " << pCD_Prey->name << " ... " << std::endl;
    for( int n = 0 ; n < parameters.ints("number_of_biofilm_cells") ; n++ ) {
      std::vector<double> position = {0,0,0};
      position[0] = Xmin + UniformRandom()*Xrange;
      position[1] = Ymin + UniformRandom()*Yrange;
      position[2] = Zmin + UniformRandom()*Zrange;
      pC = create_cell( *pCD_Prey );
      pC->assign_position( position );
    }
  }

  // place predators (probiotic_E_coli)
  Cell_Definition* pCD_Pred = find_cell_definition( "probiotic_E_coli");
  if (!pCD_Pred) {
    std::cerr << "ERROR: Could not find Cell_Definition for probiotic_E_coli during tissue setup!" << std::endl;
  } else {
    std::cout << "Placing cells of type " << pCD_Pred->name << " ... " << std::endl;
    for( int n = 0 ; n < parameters.ints("number_of_probiotic_cells") ; n++ ) {
      std::vector<double> position = {0,0,0};
      position[0] = Xmin + UniformRandom()*Xrange;
      position[1] = Ymin + UniformRandom()*Yrange;
      position[2] = Zmin + UniformRandom()*Zrange;
      pC = create_cell( *pCD_Pred );
      pC->assign_position( position );
    }
  }

  load_cells_from_pugixml();
  return;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// my_coloring_function
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<std::string> my_coloring_function( Cell* pCell )
{
  if( pCell->type == biofilm_cell_type_index ) { return { "red", "black", "darkred", "darkred" }; }
  if( pCell->type == probiotic_cell_type_index ) { return { "lime", "black", "darkgreen", "darkgreen" }; }
  return paint_by_number_cell_coloring(pCell);
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// phenotype_function (Default - likely unused)
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ) { return; }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// custom_function (Default - likely unused)
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void custom_function( Cell* pCell, Phenotype& phenotype , double dt ) { return; }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: avoid_boundaries (As provided in the file)
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void avoid_boundaries( Cell* pCell )
{
  static double Xmin = microenvironment.mesh.bounding_box[0];
  static double Ymin = microenvironment.mesh.bounding_box[1];
  static double Zmin = microenvironment.mesh.bounding_box[2];
  static double Xmax = microenvironment.mesh.bounding_box[3];
  static double Ymax = microenvironment.mesh.bounding_box[4];
  static double Zmax = microenvironment.mesh.bounding_box[5];
  static double avoid_zone = 25;
  static double avoid_speed = -0.5;
  bool near_edge = false;
  if( pCell->position[0] < Xmin + avoid_zone || pCell->position[0] > Xmax - avoid_zone ) { near_edge = true; }
  if( pCell->position[1] < Ymin + avoid_zone || pCell->position[1] > Ymax - avoid_zone ) { near_edge = true; }
  if( default_microenvironment_options.simulate_2D == false ) {
    if( pCell->position[2] < Zmin + avoid_zone || pCell->position[2] > Zmax - avoid_zone ) { near_edge = true; }
  }
  if( near_edge ) { pCell->velocity = pCell->position; pCell->velocity *= avoid_speed; }
  return;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: wrap_boundaries (Calls avoid_boundaries in the provided file)
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void wrap_boundaries( Cell* pCell ) { return avoid_boundaries( pCell ); }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: get_possible_neighbors (As provided in the file)
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<Cell*> get_possible_neighbors( Cell* pCell)
{
  std::vector<Cell*> neighbors = {};
  std::vector<Cell*>::iterator neighbor;
  std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
  for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor) { neighbors.push_back( *neighbor ); }
  std::vector<int>::iterator neighbor_voxel_index;
  std::vector<int>::iterator neighbor_voxel_index_end = pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
  for( neighbor_voxel_index = pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin(); neighbor_voxel_index!= neighbor_voxel_index_end; ++neighbor_voxel_index) {
    if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index)) continue;
    end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
    for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor) { neighbors.push_back( *neighbor ); }
  }
  return neighbors;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prey (S_aureus_biofilm) functions
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void prey_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  if (glucose_substrate_index < 0) { return; } // Safety check
  double glucose = pCell->nearest_density_vector()[glucose_substrate_index];

  // death based on glucose
  static int nNecrosis = phenotype.death.find_death_model_index( "necrosis" );
  if (nNecrosis < 0) return; // Safety check
  double necrosis_glucose_threshold = parameters.doubles("prey_necrosis_glucose_threshold");
  if( glucose < necrosis_glucose_threshold ) {
    pCell->start_death( nNecrosis );
    pCell->functions.update_phenotype = NULL;
    return;
  }

  // division based on glucose
  static Cell_Definition* pCD = find_cell_definition( "S_aureus_biofilm" );
  if (!pCD) return; // Safety check
  int cycle_start_phase_index = phenotype.cycle.current_phase_index();
  if( cycle_start_phase_index != 0 ) { return; } // Only update exit rate for G0/G1 phase (index 0)

  double base_exit_rate = pCD->phenotype.cycle.data.exit_rate(cycle_start_phase_index);
  double division_glucose_threshold = parameters.doubles("prey_division_glucose_threshold");
  double division_glucose_saturation = parameters.doubles("prey_division_glucose_saturation");
  double multiplier = 0.0;
  if (division_glucose_saturation > division_glucose_threshold) {
    multiplier = (glucose - division_glucose_threshold) / (division_glucose_saturation - division_glucose_threshold);
  }
  if (multiplier < 0.0) multiplier = 0.0;
  if (multiplier > 1.0) multiplier = 1.0;
  phenotype.cycle.data.exit_rate(cycle_start_phase_index) = base_exit_rate * multiplier;

  return;
}

void prey_custom_function( Cell* pCell, Phenotype& phenotype, double dt ) { wrap_boundaries( pCell ); }
void prey_motility_function( Cell* pCell, Phenotype& phenotype, double dt ) { return; }


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Predator (probiotic_E_coli) functions
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// --- Phenotype function modified for Conditional Neighbor Death Trigger ---
void predator_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // --- Safety Check: Ensure Indices are Valid ---
  if( glucose_substrate_index < 0 || biofilm_cell_type_index < 0 || prey_apoptosis_index < 0 ||
      dnase_substrate_index < 0 || il10_substrate_index < 0 ) // *** NEW: Check DNase/IL10 indices ***
  {
#pragma omp critical
    {
      std::cerr << "ERROR: In predator_phenotype_function (Cell ID: " << pCell->ID
        << "): Required index invalid! Skipping phenotype update." << std::endl;
      std::cerr << "       glucose_idx=" << glucose_substrate_index
        << ", biofilm_idx=" << biofilm_cell_type_index
        << ", prey_apoptosis_idx=" << prey_apoptosis_index
        << ", dnase_idx=" << dnase_substrate_index   // *** NEW ***
        << ", il10_idx=" << il10_substrate_index     // *** NEW ***
        << std::endl;
    }
    return; // Stop processing for this cell
  }

  // --- Get Parameters from XML ---
  static double glucose_threshold = parameters.doubles( "glucose_threshold");
  static double ph_threshold = parameters.doubles( "pH_threshold");
  static double current_pH = parameters.doubles( "pH_condition"); // Global pH from XML
  static double max_interaction_distance = parameters.doubles( "probiotic_interaction_distance");
  // *** NEW: Get secretion rate parameters ***
  static double dnase_secretion_rate_param = parameters.doubles( "dnase_secretion_rate" );
  static double il10_secretion_rate_param = parameters.doubles( "il10_secretion_rate" );

  // --- Check Environmental Conditions ---
  double local_glucose = pCell->nearest_density_vector()[glucose_substrate_index];
  // Condition for KILLING (requires BOTH low pH and high glucose)
  bool kill_conditions_met = ( local_glucose > glucose_threshold && current_pH < ph_threshold );
  // Condition for DNase secretion (requires low pH)
  bool dnase_condition_met = ( current_pH < ph_threshold );
  // Condition for IL10 secretion (requires high glucose)
  bool il10_condition_met = ( local_glucose > glucose_threshold );


  // *** NEW: Update Secretion Rates based on conditions ***
  // DNase secretion
  if( dnase_condition_met ) {
    phenotype.secretion.secretion_rates[dnase_substrate_index] = dnase_secretion_rate_param;
  } else {
    phenotype.secretion.secretion_rates[dnase_substrate_index] = 0.0;
  }
  // IL10 secretion
  if( il10_condition_met ) {
    phenotype.secretion.secretion_rates[il10_substrate_index] = il10_secretion_rate_param;
  } else {
    phenotype.secretion.secretion_rates[il10_substrate_index] = 0.0;
  }


  // --- Conditional Behavior: Trigger Apoptosis in Nearby Prey (if kill conditions met) ---
  if( kill_conditions_met ) // Killing still requires BOTH conditions
  {
    std::vector<Cell*> neighbors = get_possible_neighbors( pCell);
    for( Cell* pNeighbor : neighbors ) {
      if (!pNeighbor || pNeighbor == pCell) { continue; }

      if( pNeighbor->type == biofilm_cell_type_index && pNeighbor->phenotype.death.dead == false ) {
        std::vector<double> displacement = pNeighbor->position; displacement -= pCell->position;
        double distance_squared = norm_squared( displacement );
        double interaction_radius = pCell->phenotype.geometry.radius + pNeighbor->phenotype.geometry.radius + max_interaction_distance;
        double interaction_radius_squared = interaction_radius * interaction_radius;

        if( distance_squared < interaction_radius_squared ) {
          pNeighbor->start_death( prey_apoptosis_index );
          // Optional logging (use omp critical if multi-threaded)
          // #pragma omp critical
          // { std::cout << "Time: " << PhysiCell_globals.current_time << " Probiotic " << pCell->ID << " triggered apoptosis in Biofilm " << pNeighbor->ID << std::endl; }
          break; // Kill only one neighbor per predator per time step
        }
      }
    }
  }
  // --- End Conditional Killing Logic ---


  // --- Optional: Energy / Death logic ---
  /* ... (remains commented out unless needed) ... */

  return; // End of function
}
void predator_phenotype_function1( Cell* pCell, Phenotype& phenotype, double dt )
{
  // --- Safety Check: Ensure Indices are Valid ---
  if( glucose_substrate_index < 0 || biofilm_cell_type_index < 0 || prey_apoptosis_index < 0 ) // Check apoptosis index too
  {
#pragma omp critical
    {
      std::cerr << "ERROR: In predator_phenotype_function (Cell ID: " << pCell->ID
        << "): Required index invalid! Skipping phenotype update." << std::endl;
      std::cerr << "       glucose_idx=" << glucose_substrate_index
        << ", biofilm_idx=" << biofilm_cell_type_index
        << ", prey_apoptosis_idx=" << prey_apoptosis_index << std::endl;
    }
    return; // Stop processing for this cell
  }

  // --- Get Parameters from XML ---
  // Use parameters.doubles("name") syntax first, matching parameters.ints("name") usage
  static double glucose_threshold = parameters.doubles( "glucose_threshold");
  static double ph_threshold = parameters.doubles( "pH_threshold");
  static double current_pH = parameters.doubles( "pH_condition"); // Global pH from XML
  static double max_interaction_distance = parameters.doubles( "probiotic_interaction_distance");

  // --- Check Environmental Conditions ---
  double local_glucose = pCell->nearest_density_vector()[glucose_substrate_index];
  bool conditions_met = ( local_glucose > glucose_threshold && current_pH < ph_threshold );

  // --- Conditional Behavior: Trigger Apoptosis in Nearby Prey ---
  if( conditions_met )
  {
    // Conditions met: Find nearby prey cells and trigger their apoptosis

    // Use the existing neighbor function from the provided file
    std::vector<Cell*> neighbors = get_possible_neighbors( pCell);

    for( Cell* pNeighbor : neighbors ) // Iterate using range-based for loop
    {
      // Safety check & ignore self
      if (!pNeighbor || pNeighbor == pCell) { continue; }

      // Check if the neighbor is a biofilm cell AND is currently alive
      if( pNeighbor->type == biofilm_cell_type_index && pNeighbor->phenotype.death.dead == false )
      {
        // Check distance
        std::vector<double> displacement = pNeighbor->position;
        displacement -= pCell->position;
        double distance_squared = norm_squared( displacement ); // Use squared distance

        // Interaction radius = sum of radii + max_interaction_distance
        double interaction_radius = pCell->phenotype.geometry.radius + pNeighbor->phenotype.geometry.radius + max_interaction_distance;
        double interaction_radius_squared = interaction_radius * interaction_radius;

        // If within interaction radius, trigger apoptosis
        if( distance_squared < interaction_radius_squared )
        {
          // Trigger apoptosis in the neighbor cell
          pNeighbor->start_death( prey_apoptosis_index );

          // Optional: Logging
#pragma omp critical
          {
            std::cout << "Time: " << PhysiCell_globals.current_time << " Probiotic " << pCell->ID
              << " triggered apoptosis in Biofilm " << pNeighbor->ID << std::endl;
          }

          // Break the inner loop after triggering death for one neighbor
          // This prevents one predator from killing multiple prey in the exact same time step.
          // Adjust this if different behavior is desired.
          break;

        } // end if within distance
      } // end if neighbor is correct type and alive
    } // end for neighbors loop
  } // end if conditions_met
    // else: Conditions are not met, so this probiotic does not actively kill neighbors in this step.

    // --- Optional: Original Energy / Death logic (Can be kept or removed) ---
  /*
     static double decay_rate = parameters.doubles("predator_energy_decay_rate");
     if( pCell->custom_data.find("energy") != pCell->custom_data.end() ) {
     pCell->custom_data["energy"] /= (1.0 + dt*decay_rate);
     } else { pCell->custom_data["energy"] = parameters.doubles("predator_initial_energy"); }

     static int nNecrosis = phenotype.death.find_death_model_index( "necrosis" );
     if( nNecrosis >= 0 && pCell->custom_data["energy"] < parameters.doubles("predator_energy_death_threshold") )
     { pCell->start_death( nNecrosis ); pCell->functions.update_phenotype = NULL; return; }
     */

  return; // End of function
}


void predator_custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  wrap_boundaries( pCell );
  return;
}

void predator_motility_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  if (glucose_substrate_index >= 0 && parameters.bools("predator_chemotaxis_enabled") ) // Check if enabled
  {
    phenotype.motility.migration_bias_direction = pCell->nearest_gradient(glucose_substrate_index);
    normalize( &(phenotype.motility.migration_bias_direction) );
  }
  else { phenotype.motility.migration_bias_direction = {0,0,0}; }
  return;
}
