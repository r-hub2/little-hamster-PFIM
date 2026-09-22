# Stable accessors for Optimization / Design tests (avoid raw property names).

.opt_outputs = function( opt ) {
  PFIM:::.pfimAlgoOutputs( opt )
}

.opt_algo = function( opt ) {
  PFIM:::.getOptimizationAlgorithm( opt )
}

.opt_eval = function( opt ) {
  PFIM:::.getOptimalEval( opt )
}

.opt_init_eval = function( opt ) {
  PFIM:::.getInitialEval( opt )
}

.opt_mult_out = function( opt ) {
  prop( .opt_algo( opt ), "multiplicativeAlgorithmOutputs" )
}

.opt_fw_out = function( opt ) {
  prop( .opt_algo( opt ), "FedorovWynnAlgorithmOutputs" )
}

.opt_algo_status = function( opt ) {
  prop( .opt_algo( opt ), "optimizerOutputs" )$algorithmOutput
}

.project_design = function( project, i = 1L ) {
  projectProp( project, "designs" )[[ i ]]
}

.design_arms = function( design ) {
  prop( design, "arms" )
}

.first_arm = function( project ) {
  .design_arms( .project_design( project ) )[[ 1L ]]
}

.eval_design = function( evaluation, i = 1L ) {
  prop( evaluation, "designs" )[[ i ]]
}

.fim_matrix = function( fim ) {
  prop( fim, "fisherMatrix" )
}
