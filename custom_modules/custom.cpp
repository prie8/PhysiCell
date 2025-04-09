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
int dnase_substrate_index = -1;
int il10_substrate_index = -1;
int biofilm_cell_type_index = -1;
int probiotic_cell_type_index = -1;
int prey_apoptosis_index = -1; // Index for the prey's apoptosis death model

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// create_cell_types function - Uses parenthesis parameter access
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_cell_types( void )
{
  // set the random seed using parenthesis notation
  SeedRandom( parameters.ints("random_seed") ); // Read from user_parameters

  initialize_default_cell_definition();
  cell_defaults.functions.volume_update_function = standard_volume_update_function;
  cell_defaults.functions.update_velocity = standard_update_cell_velocity;
  cell_defaults.functions.update_phenotype = NULL; // Phenotype functions defined per type
  cell_defaults.functions.custom_cell_rule = NULL; // Custom rules defined per type
  cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
  cell_defaults.functions.calculate_distance_to_membrane = NULL;
  // Assign the custom coloring function to the default cell definition
  // This ensures it's used for all cell types unless overridden specifically
  // cell_defaults.functions.coloring_function = my_coloring_function; // Assign coloring function

  initialize_cell_definitions_from_pugixml();

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

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// setup_microenvironment function - Unchanged
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
  glucose_substrate_index = microenvironment.find_density_index("glucose");
  dnase_substrate_index = microenvironment.find_density_index("dnase");
  il10_substrate_index = microenvironment.find_density_index("il10");
  if (glucose_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'glucose'!" << std::endl; }
  if (dnase_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'dnase'!" << std::endl; }
  if (il10_substrate_index < 0) { std::cerr << "ERROR: Could not find microenvironment density 'il10'!" << std::endl; }
  return;
}


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// setup_tissue function - Uses parenthesis parameter access
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void setup_tissue( void )
{
  Cell* pC;
  // Define radii for placement - Read parameters using parenthesis notation
  double wound_radius = parameters.doubles( "wound_radius" );
  double biofilm_placement_radius = parameters.doubles( "initial_biofilm_radius" );
  double predator_min_placement_radius = biofilm_placement_radius + parameters.doubles( "predator_placement_buffer" );
  double predator_max_placement_radius = wound_radius - parameters.doubles( "predator_boundary_buffer" );

  // place prey (S_aureus_biofilm) in a central disc
  Cell_Definition* pCD_Prey = find_cell_definition( "S_aureus_biofilm");
  if (!pCD_Prey) { std::cerr << "ERROR: Could not find Cell_Definition for S_aureus_biofilm during tissue setup!" << std::endl; }
  else {
    std::cout << "Placing biofilm cells in central disc (radius " << biofilm_placement_radius << ") ... " << std::endl;
    int num_biofilm = parameters.ints("number_of_biofilm_cells");
    for( int n = 0 ; n < num_biofilm ; n++ ) {
      std::vector<double> position = {0,0,0};
      double r = sqrt(UniformRandom()) * biofilm_placement_radius; double theta = UniformRandom() * 2.0 * M_PI;
      position[0] = r * cos(theta); position[1] = r * sin(theta);
      pC = create_cell( *pCD_Prey ); pC->assign_position( position );
    }
  }
  // place predators (probiotic_E_coli) outside the biofilm disc, within wound radius
  Cell_Definition* pCD_Pred = find_cell_definition( "probiotic_E_coli");
  if (!pCD_Pred) { std::cerr << "ERROR: Could not find Cell_Definition for probiotic_E_coli during tissue setup!" << std::endl; }
  else {
    if (predator_max_placement_radius <= predator_min_placement_radius) {
      predator_max_placement_radius = predator_min_placement_radius + 10.0;
      std::cout << "Warning: Predator max placement radius adjusted to " << predator_max_placement_radius << std::endl;
    }
    std::cout << "Placing probiotic cells outside disc (radius " << predator_min_placement_radius << " to " << predator_max_placement_radius << ") ... " << std::endl;
    int num_probiotic = parameters.ints("number_of_probiotic_cells");
    for( int n = 0 ; n < num_probiotic ; n++ ) {
      std::vector<double> position = {0,0,0};
      double r_squared = UniformRandom() * (predator_max_placement_radius*predator_max_placement_radius - predator_min_placement_radius*predator_min_placement_radius) + predator_min_placement_radius*predator_min_placement_radius;
      double r = sqrt(r_squared); double theta = UniformRandom() * 2.0 * M_PI;
      position[0] = r * cos(theta); position[1] = r * sin(theta);
      pC = create_cell( *pCD_Pred ); pC->assign_position( position );
    }
  }
  load_cells_from_pugixml();
  return;
} // End of setup_tissue


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// my_coloring_function - Assumes DNase=true, reads IL10 state from custom data
// Uses parenthesis parameter access for pH check
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<std::string> my_coloring_function( Cell* pCell )
{
  std::vector<std::string> output = paint_by_number_cell_coloring(pCell);

  // Biofilm is Red
  if( pCell->type == biofilm_cell_type_index )
  { output[0] = "red"; output[1] = "black"; output[2] = "darkred"; output[3] = "black"; }
  // Probiotic color
  else if( pCell->type == probiotic_cell_type_index )
  {
    // --- Determine Secretion States ---
    // Read IL10 state from custom data (this seemed reliable)
    bool secreting_il10 = pCell->custom_data["secreting_il10"] > 0.5; // Assumes key exists after first phenotype update

    // Assume DNase condition is always TRUE based on constant pH < threshold
    // If parameter reading becomes reliable inside this function, this could be re-enabled:
    double ph_threshold = parameters.doubles( "pH_threshold");
    double current_pH = parameters.doubles( "pH_condition");
    bool secreting_dnase = ( current_pH < ph_threshold );
    // bool secreting_dnase = true;
    output[0] = "lime"; output[1] = "black"; output[2] = "darkgreen"; output[3] = "black";

    // Color based on assumed DNase state and stored IL10 state
    if( secreting_dnase && secreting_il10 ) { output[0] = "lime"; }
    else if ( secreting_dnase ) { output[0] = "yellow"; }
    else if ( secreting_il10 ) { output[0] = "cyan"; }
    else { output[0] = "magenta"; }

    output[2] = output[0]; // Match nucleus color
  }
  return output;
}


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// phenotype_function (Default - unused) - Unchanged
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ) { return; }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// custom_function (Default - unused) - Unchanged
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void custom_function( Cell* pCell, Phenotype& phenotype , double dt ) { return; }


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: Confine cells to a circular wound area - Uses parenthesis parameter access
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void confine_to_wound_area( Cell* pCell, double wound_radius ) // Pass radius
{
  if( wound_radius <= 0.0 ) { return; }
  double dist_from_center_squared = pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1];
  double wound_radius_squared = wound_radius * wound_radius;

  if( dist_from_center_squared > wound_radius_squared ) {
    std::vector<double> push_direction = { -pCell->position[0], -pCell->position[1], 0.0 };
    double norm_factor = norm(push_direction);
    if (norm_factor > 1e-16) {
      push_direction[0] /= norm_factor; push_direction[1] /= norm_factor;
      double push_speed = pCell->phenotype.motility.migration_speed;
      if (push_speed < 0.1) { push_speed = 0.5; } // Ensure minimum push speed
      pCell->velocity[0] = push_speed * push_direction[0];
      pCell->velocity[1] = push_speed * push_direction[1];
      if( default_microenvironment_options.simulate_2D == false ) { pCell->velocity[2] = 0; }
    } else { pCell->velocity = {0,0,0}; }
  }
}


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: avoid_boundaries - Unchanged
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void avoid_boundaries( Cell* pCell ) { /* ... unchanged ... */ }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: wrap_boundaries - Unchanged
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void wrap_boundaries( Cell* pCell ) { return avoid_boundaries( pCell ); }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Helper: get_possible_neighbors - REMOVED as unused
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// std::vector<Cell*> get_possible_neighbors( Cell* pCell) { /* ... removed ... */ }


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prey (S_aureus_biofilm) functions - Uses parenthesis parameter access
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void prey_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  if (glucose_substrate_index < 0) { return; }
  double glucose = pCell->nearest_density_vector()[glucose_substrate_index];

  // death based on glucose
  static int nNecrosis = phenotype.death.find_death_model_index( "necrosis" );
  if (nNecrosis < 0) return;
  // Use parenthesis notation for parameter access - REMOVED STATIC
  double necrosis_glucose_threshold = parameters.doubles("prey_necrosis_glucose_threshold");
  if( glucose < necrosis_glucose_threshold ) {
    pCell->start_death( nNecrosis );
    pCell->functions.update_phenotype = NULL;
    return;
  }

  // division based on glucose
  static Cell_Definition* pCD = find_cell_definition( "S_aureus_biofilm" );
  if (!pCD) return;
  int cycle_start_phase_index = phenotype.cycle.current_phase_index();
  if( cycle_start_phase_index != 0 ) { return; }

  double base_exit_rate = pCD->phenotype.cycle.data.exit_rate(cycle_start_phase_index);
  // Use parenthesis notation for parameter access - REMOVED STATIC
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

// MODIFIED custom function - Uses parenthesis parameter access
void prey_custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // Read parameter using parenthesis notation - REMOVED STATIC
  double wound_radius = parameters.doubles("wound_radius");
  confine_to_wound_area( pCell, wound_radius );
}
void prey_motility_function( Cell* pCell, Phenotype& phenotype, double dt ) { return; } // No changes needed


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Predator (probiotic_E_coli) functions - Uses parenthesis parameter access
// Sets only IL10 custom data, uses built-in neighbor search
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void predator_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // --- Safety Check ---
  if( glucose_substrate_index < 0 || biofilm_cell_type_index < 0 || prey_apoptosis_index < 0 ||
      dnase_substrate_index < 0 || il10_substrate_index < 0 )
  { /* ... error ... */ return; }

  // --- Get Parameters --- using parenthesis syntax, removed static
  double glucose_threshold = parameters.doubles( "glucose_threshold");
  double ph_threshold = parameters.doubles( "pH_threshold");
  double current_pH = parameters.doubles( "pH_condition");
  double max_interaction_distance = parameters.doubles( "probiotic_interaction_distance");
  double dnase_secretion_rate_param = parameters.doubles( "dnase_secretion_rate" );
  double il10_secretion_rate_param = parameters.doubles( "il10_secretion_rate" );
  double kill_probability = parameters.doubles( "probiotic_kill_probability" );

  // --- Check Conditions ---
  double local_glucose = pCell->nearest_density_vector()[glucose_substrate_index];
  bool dnase_condition_met = ( current_pH < ph_threshold );
  bool il10_condition_met = ( local_glucose > glucose_threshold );
  bool kill_conditions_met = dnase_condition_met;

  // --- Update Secretion Rates ---
  if (dnase_substrate_index < phenotype.secretion.secretion_rates.size()) {
    phenotype.secretion.secretion_rates[dnase_substrate_index] = dnase_condition_met ? dnase_secretion_rate_param : 0.0;
  }
  if (il10_substrate_index < phenotype.secretion.secretion_rates.size()) {
    phenotype.secretion.secretion_rates[il10_substrate_index] = il10_condition_met ? il10_secretion_rate_param : 0.0;
  }

  // --- Conditional Killing --- using built-in neighbor search
  if( kill_conditions_met ) {
    // Use PhysiCell's standard function (store by value)
    std::vector<Cell*> neighbors = pCell->nearby_interacting_cells(); // USE BUILT-IN
    for( Cell* pNeighbor : neighbors ) {
      if (!pNeighbor || pNeighbor == pCell) { continue; } // Skip self/null
      if( pNeighbor->type == biofilm_cell_type_index && pNeighbor->phenotype.death.dead == false ) { // Check type and liveness
        std::vector<double> displacement = pNeighbor->position; displacement -= pCell->position;
        double distance_squared = norm_squared( displacement );
        double interaction_radius = pCell->phenotype.geometry.radius + pNeighbor->phenotype.geometry.radius + max_interaction_distance;
        double interaction_radius_squared = interaction_radius * interaction_radius;
        if( distance_squared < interaction_radius_squared ) { // Check distance
// #pragma omp critical
//           { std::cout << "KILL TRIGGERED..." << std::endl; }
           if( UniformRandom() < kill_probability ) {
          pNeighbor->start_death( prey_apoptosis_index ); // Trigger death
                                                          // Optional: Add kill log here if needed
          break; // Kill only one neighbor per step
        }
        }
      }
    }
  }
  // --- End Conditional Killing Logic ---

  // --- Store *only IL10* secretion state in custom data ---
  pCell->custom_data["secreting_il10"] = il10_condition_met ? 1.0 : 0.0;

  return;
}

// MODIFIED custom function - Uses parenthesis parameter access
void predator_custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // Read parameter using parenthesis notation - removed static
  double wound_radius = parameters.doubles("wound_radius");
  confine_to_wound_area( pCell, wound_radius );
}

// MODIFIED motility function - Uses parenthesis parameter access, stronger confinement
void predator_motility_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // Get parameters for vicinity confinement (using parenthesis syntax) - removed static
  double initial_biofilm_radius = parameters.doubles("initial_biofilm_radius");
  double vicinity_buffer = parameters.doubles("probiotic_max_distance_buffer");
  double confinement_speed = parameters.doubles("probiotic_confinement_speed");
  double max_allowed_radius = initial_biofilm_radius + vicinity_buffer;
  bool chemotaxis_enabled = parameters.bools("predator_chemotaxis_enabled");

  double dist_from_center_squared = pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1];
  double max_allowed_radius_squared = max_allowed_radius * max_allowed_radius;

  if( dist_from_center_squared > max_allowed_radius_squared && max_allowed_radius > 0 )
  {
    // Confine: Force movement back towards center by setting velocity directly
    phenotype.motility.is_motile = true;
    std::vector<double> direction_to_center = { -pCell->position[0], -pCell->position[1], 0.0 };
    normalize( direction_to_center );
    pCell->velocity = confinement_speed * direction_to_center; // Use confinement speed
    phenotype.motility.migration_bias_direction = {0,0,0}; // Clear bias direction as velocity is set
  }
  else
  {
    // Normal Behavior: Chemotaxis towards glucose (within vicinity)
    phenotype.motility.is_motile = true;
    // Restore base bias value from phenotype? Not needed if only direction is set.
    // phenotype.motility.migration_bias = pCell->phenotype.motility.migration_bias;

    if (glucose_substrate_index >= 0 && chemotaxis_enabled )
    {
      phenotype.motility.migration_bias_direction = pCell->nearest_gradient(glucose_substrate_index);
      normalize( &(phenotype.motility.migration_bias_direction) );
      // Let PhysiCell core use the XML migration_bias value with this direction
    }
    else
    {
      phenotype.motility.migration_bias_direction = {0,0,0}; // No bias direction
    }
  }
  return;
}
