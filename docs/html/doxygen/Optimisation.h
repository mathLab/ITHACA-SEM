/**
 * @page Optimisation Optimisation
 * <h2>Global Optimisation</h2>
 * It is possible to specify the following 2 optimisation flags for different
 * operators:    
 * - DO_GLOBAL_MAT_OP
 * If TRUE (VALUE = 1), the globally assembled system matrix will be used to
 * evaluate the operator. If FALSE (VALUE = 0), the operator will be evaluated
 * elementally.
 * - DO_BLOCK_MAT_OP
 * If TRUE (VALUE = 1), the elemental evaluation will be done using the
 * elemental/local matrices (which are all concatenated in a block matrix,
 * hence the name). If FALSE (VALUE = 0), the elemental evaluation will be done
 * using the sum-factorisation technique.
 *
 * The optimal configuration of these parameters depends on
 * - the computer environment (i.e. processor, BLAS library, ...)
 * - the discretisation (the mesh-size but mainly the POLYNOMIAL ORDER)
 *
 * <h3>Choosing parameters</h3>
 * To set these parameters use the following rules of thumb:
 * - The global matrix approach is the most efficient option
 *     only for very low expansion order (P=1,P=2). Never set 
 *     DO_GLOBAL_MAT_OP to true for an expansion order P>4 as it 
 *     quickly becomes very expensive and you might run
 *     out of memory.
 * - The most efficient way of elementally evaluating an
 *     operator depends on the complexity of the operator
 *     in the sum-factorisation approach.
 *     -# For simple operators (such as BwdTrans and IProductWRTBase)
 *       - set DO_BLOCK_MAT_OP to TRUE for low orders (e.g. P<=4)
 *       - set DO_BLOCK_MAT_OP to FALSE in the other case
 *     -# For more complex operators (such as MassMatrixOp)
 *       - set DO_BLOCK_MAT_OP to TRUE for low and intermediate orders (e.g. P<=8)
 *       - set DO_BLOCK_MAT_OP to FALSE in the other case
 *     -# For very complex operators (such as HelmholtzMatrixOp)
 *       always set DO_BLOCK_MAT_OP to TRUE
 *     In general, the break-even point between the elemental matrix approach 
 *     and the sum-factorisation technique is higher for triangular meshes than
 *     for quadrilateral meshes. That is why you for exmaple may want to set 
 *     following flags for the IProductWRTBase operator:
 *       (quadrilateral mesh) set DO_BLOCK_MAT_OP to TRUE if P<=4
 *       (triangular mesh)    set DO_BLOCK_MAT_OP to TRUE if P<=6
 *
 *   Also note that, as opposed to the ELEMENTAL optimisation parameters,
 *   the GLOBAL parameters are problem dependent. In order to optimally set
 *   these parameters below, ideally, a 'self-tuning' optimisation suite
 *   which performs a series of test-runs based upon the mesh file above
 */
