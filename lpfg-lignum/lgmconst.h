/* ***************************************************** */
/*                  LPFG-LIGNUM constants                */
/* ***************************************************** */

#ifndef RNDSEED
/* *** Random number generator seed [int]
   0 stands for automatic and all other values are actual fixed seeds.
*/
#define RNDSEED 1
#endif

#ifndef VERB
/* *** Verbose level [int]
   1 is for "Year: ..., M-P = ..., Lambda = ... " output.
   0 is for turning off the messages, however, warnings and errors will be there
   anyway.
*/
#define VERB 1
#endif

#ifndef YEARS
/* *** Number of years to simulate [int]
   The time step of LIGNUM is 1 year. This eventually means the number of
   iterations.
*/
#define YEARS 10
#endif

#ifndef AF
/* *** Foliage mass coeff [kg/m^2]
   A factor to the segments surface area, determining the foliage mass of the
   segment. Namely, foliage mass = AF * 2 * PI * Len * Rad.
 */
#define AF 1.3
#endif

#ifndef AR
/* *** Root increment coeff []
   The new foliage induces the root growth.
   Root_mass_increment = AR * fol_mass_increment. 
 */
#define AR 0.5
#endif

#ifndef SR
/* *** Root senescence rate [year^-1]
   SR portion of the root mass dies out annually
 */
#define SR 0.33
#endif

#ifndef SS
/* *** Sapwood senescence [year^-1]
   SS portion of the sapwood cross sectional area turns into heartwood annually. The
   variation: 0.2 - 1.0 [Perttunen 1998], fractional value, i.e. 0 = 0%, 1.0 = 100%.
 */
#define SS 0.2
#endif

#ifndef MF
/* *** Foliage respiration rate []
   Portion of the foliage mass that is converted to the mass-energy consumed for the
   respiration purposes. Variation: up to 1.0 (fractional value).
 */
#define MF 0.2
#endif

#ifndef MS
/* *** Sapwood respiration rate []
   Portion of the sapwood mass that is converted to the mass-energy consumed for the
   respiration purposes. Variation: up to 1.0 (fractional value).
 */
#define MS 0.024
#endif

#ifndef MR
/* *** Roots respiration rate []
   Portion of the root mass that is converted to the mass-energy consumed for the
   respiration purposes. Variation: up to 1.0 (fractional value).
 */
#define MR 0.24
#endif

#ifndef KSI
/* *** Initial heartwood []
   The initial heartwood propoprtion of a segment calculated as its cross sectional
   area radius, i.e. R_h = sqrt(KSI) * R_segment. Variation: 0.6-0.7 [Sievanen 2008]
 */
#define KSI 0.6
#endif

#ifndef RHO
/* *** Wood density [kg/m^3]
   Used to calculate masses. Variation: 200-400 [Perttunen 1998, Sievanen 2008].
 */
#define RHO 400.0
#endif

#ifndef NEEDLELEN
/* *** Needle length [m]
   Used only for demonstration purposes in LPFG-LIGNUM.
 */
#define NEEDLELEN 0.03
#endif

#ifndef NEEDLEANG
/* *** Needle inclination angle [rad]
   Used only for demonstration purposes in LPFG-LIGNUM.
 */
#define NEEDLEANG (M_PI/4)
#endif

#ifndef NUMNEEDLECIRC
/* *** Number of needles in the circumference [int]
   Used only for demonstration purposes in LPFG-LIGNUM.
 */
#define NUMNEEDLECIRC 5
#endif

#ifndef NUMNEEDLELONG
/* *** Number of needles along a segment [int]
   Used only for demonstration purposes in LPFG-LIGNUM.
 */
#define NUMNEEDLELONG 10
#endif

#ifndef NEEDLEDENMIN
/* *** Minimum distance between needles [m]
   Determines the density of needles on a segment. Used only for demonstration
   purposes in LPFG-LIGNUM.
 */
#define NEEDLEDENMIN 0.02
#endif

#ifndef NEEDLEDENMAX
/* *** Maximum distance between needles [m]
   Determines the density of needles on a segment. Used only for demonstration
   purposes in LPFG-LIGNUM.
 */
#define NEEDLEDENMAX 0.05
#endif

#ifndef NEEDLEDIAM
/* *** Diameter of a needle [m]
   Used only for demonstration purposes in LPFG-LIGNUM.
 */
#define NEEDLEDIAM 0.002
#endif

#ifndef RF
/* *** Foliage induced radius [m]
   Foliage creates an addition to the segment's radius. This is used in spatial
   arrangements of the segments.
 */
#define RF (NEEDLELEN*cos(M_PI/2 - NEEDLEANG))
#endif

#ifndef BRAANG
/* *** Branching angle [deg]
   Initial inclination angle of the child segment (offset). Variation: normally,
   around 45 deg, but can be anything ~30 - ~60, depending on the final structure of
   the tree.
 */
#define BRAANG 45
#endif

#ifndef ANGINCR
/* *** Branching angle increment [deg]
   When new segment added to the tip of a branch, the branch gets bent. That is, all
   its segments get additional inclination to the parent segment, the branch emanates
   from. Variation: normally around 5-15 deg, but can be anything positive, depends
   highly on BRAANG and MAXANG.
 */
#define ANGINCR 10.0
#endif

#ifndef MAXANG
/* *** Maximum branching angle [deg]
   While bending the segments in a branch cannot get the inclination angle larger
   than this value. Variation: normally around 90 deg, can be slightly less or more
   than this, e.g. ~80 - ~105 deg.
 */
#define MAXANG 95.0
#endif

#ifndef INITROTANG
/* *** Rotation angle [deg]
   The first lateral buds (giving branching segment) are rotated ralative to each
   other by this angle around the parent segment in subsequent branching points
   (144 deg in Perttunen et al. 1998). The 144 deg value results in every 6th year
   first lateral buds emanating from the same axis point to the same direction. Not
   really crucial.
 */
#define INITROTANG 144.0
#endif

#ifndef FLA
/* *** Length function coeff []
   Determines at which photosynthetic ratio (i_p) the function becomes zero. It is a
   linear function determining the relative segment length depending on the relative
   insolation conditions i_p for the segment. It is 1.0 at i_p = 1.0 and 0.0 at i_p =
   FLA. Variation: 0.2-0.4 [Perttunen 1998, Sievanen 2008].
 */
#define FLA 0.4
#endif

#ifndef FLBUDS
/* *** Length function threshold for lateral buds []
   The minimum value of the length function that causes a segment to produce lateral
   buds according to the pre-defined empirial function. If the length function value
   is less than FLBUDS, no lateral buds are produced. Variation: 0.0 - 0.6 [Sievanen
   2008]. 
 */
#define FLBUDS 0.0
#endif

#ifndef Q
/* *** Apical dominance parameter []
   The increase in the order (w) of segments causes them to shorten corresponding to
   Q. Namely, the new segment length = Lambda * (1 - w*Q) * length_function (Lambda
   is a iteratively found parameter, from the balance equation). Variation: 0.0-1.0,
   1.0 is a full apical dominance (no lateral growth) whereas 0.0 is no apical
   dominance (bushy tree with no distinguished stem).
 */
#define Q 0.2
#endif

#ifndef MINLEN
/* *** Minimum length allowed for a segment [m]
   The new segments shorter than this are not grown. Variation: 0.0 - 0.02 [Sievanen
   2008]. 
 */
#define MINLEN 0.0
#endif

#ifndef PR
/* *** Photosynthetic efficiency [kg/mJ]
   Amount of the photosynthates (kg) per mJ of the intercepted radiation. Determines
   how efficient is the photosynthesis. Variation: 0.0006 - 0.001 [Perttunen 1998,
   Sievanen 2008].
 */
#define PR 0.001
#endif

#ifndef LR
/* *** Radius to length ratio []
   The fixed ration between the radius and length of the newly grown
   segment. Variation: 0.008 - 0.012 [Sievanen 2008].
 */
#define LR 0.01
#endif

#ifndef LRMIN
/* *** Minimum LR []
   Minimum LR ratio for the worst insolation conditions. (Mentioned only in Sievanen
   et al. 2008, NOT used here at the moment). The linear interpolation should be used
   between LRMIN and LRMAX.
 */
#define LRMIN 0.008
#endif

#ifndef LRMAX
/* *** Maximum LR []
   Maximum LR ratio for the best insolation conditions. (Mentioned only in Sievanen
   et al. 2008, NOT used here at the moment). The linear interpolation should be used
   between LRMIN and LRMAX.
 */
#define LRMAX 0.012
#endif

#ifndef INITLEN
/* *** Length of the very first segment [m]
   The seed segment (first year) has the pre-defined length. Variation: anything
   close to 0.3-0.5.
 */
#define INITLEN 0.5
#endif

#ifndef FOLTHR
/* *** NO-foliage threshold [kg]
   The value, which is considered to be practically zero for the foliage mass of a
   segment. Variation: must be quite small.
 */
#define FOLTHR 1e-08
#endif

#ifndef SHEDYEAR
/* *** The youngest segment age in a branch to be shed [int]
   The age criterium for the youngest segment in a branch to be shed. That is, the
   tip segment of the branch must have at least SHEDYEAR age for the whole branch to
   be shed. Note, that the branch will not be shed if it still has any foliage left,
   disregarding this value. Variation: anything <= 5 will not have an effect, since
   the default [Perttunen 1998, Sievanen 2008] empirical shedding functions remove
   completely foliage at the 5th year.
 */
#define SHEDYEAR 0
#endif

#ifndef SHEDDIST
/* *** Stochastic shedding distribution type [int]
   The shedding of the segments in a 1st order (only) branch can be stochastic if
   this value is not zero. That is, the number of segments to shed starting from
   the tip of the branch is determined randomly and the rest is left in the final
   structure of the tree. Values: 0 - shed everything, 1 - uniformly determined
   number of segments to shed, 2 - poisson determined.
*/
#define SHEDDIST 0
#endif

#ifndef SHEDMU
/* *** Mean for the Poisson shedding [float] <= 10
   In case of Poisson shedding, the mean for the Poisson distribution. Note the mean
   cannot be larger than 10. If the shedding is due to the stochastic reasons (wind,
   animals, snow fall etc.) the Poisson shedding is a good approximation to this,
   since the tip segments get shed easier, than the base ones.
 */
#define SHEDMU 2.0
#endif

#ifndef BRCRITLR
/* *** Critical ratio of the branche length to base radius []
   To imitate effect of the stochastic shedding and environmental shading to suppres
   lateral growth. NOT used in the current version.
 */
#define BRCRITLR 50000
#endif

#ifndef CRITRUPT
/* *** Critical value for the rupture modulus [kg/m^2]
   The rupture modulus describe the stability of the wood to the rupture stress.
   This value is taken from databases for the live pine tree branches and seems
   to be very high for any branch to break and shed in the normal conditions.
   NOT used in the current version
 */
#define CRITRUPT 83300// [N/m^2] OR crit MOR=(12080*698),[kg/m^2]
#endif

#ifndef PIPECOEFF
/* *** Pipe model requirement coefficient []
   The coefficient for the pipe model determining how much of the sapwood area of a
   lateral (child) segment contributes to that of the parent segment. Softening of
   the pipe model requirements. It is a single value for all orders.
   Variation: 0.75 - 1.0 [Sievanen 2008].
   */
#define PIPECOEFF 1.0
#endif

#ifndef SSSTART
/* *** Sapwood senescence start [int]
   The year (number of iteration) after which the sapwood senescence
   starts. Variation: can be 15 [Sievanen 2008].
*/
#define SSSTART 0
#endif

#ifndef GAMMASD
/* *** Horizontal curvature [deg]
   The horizontal contribution to the curvature of the branches in space. The
   XY-plane angle GAMMA is following the normal distribution
   N(0.0,GAMMASD). Variation: not really big, since produces unrealistic curvature,
   i.e. > 90 will turn the branches back (quite non-realistic), 5.0-10.0 deg maintain
   the desired level of reality.
*/
#define GAMMASD 5.0
#endif

#ifndef ZETASD
/* *** Vertical curvature [deg]
   The vertical contribution to the curvature of the branches in space. The
   XZ-plane angle ZETA is following the normal distribution
   N(ANGINCR,ZETASD). Variation: 5.0-10.0 is good enough.
*/
#define ZETASD 5.0
#endif

#ifndef ENVSHRAD
/* *** Environmental shadow [m]
   The shadow is also exerted by the neighboring trees (environment). This shadow is
   mimicked simply by not allowing the lateral branches grow beyong the circular
   boundary identified by the center at the tree seed segment and the radius
   ENVSHRAD. Variation: depending on the forest density, 1.0 - 4.0 is quite
   realistic. 
*/
#define ENVSHRAD 1.0
#endif

#ifndef GSA
/* *** Size of the voxel grid/box [int]
   Size of the voxel grid for each of the directions, for the individual directions
   see the parameters GSX, GSY, and GSZ below. Variation: depending on the VSA must
   cover whole tree in space, otherwise, the fatal error and simulation is wasted.
*/
#define GSA 600
#endif

#ifndef GSX
/* *** Size of the voxel box along X-axis [int]
 */
#define GSX GSA
#endif

#ifndef GSY
/* Size of the voxel box along Y-axis (up direction in LPFG) [int]
 */
#define GSY GSA
#endif

#ifndef GSZ
/* Size of the voxel box along Z-axis [int]*/
#define GSZ GSA
#endif

#ifndef VSA
/* *** Physical size of a voxel [m]
   The physical size of a voxel in the voxel grid for each of the directions. For the
   individual direction see the parameter VSX, VSY, and VSZ below. The shadow
   propagation model requires VSX = VSY = VSZ for better accuracy. Variation: 0.01 -
   0.05 is quite good choice, larger values introduce unnecessary discreteness to the
   space, smaller ones intensify the computation.
 */
#define VSA 0.01
#endif

#ifndef VSX
/* *** Physical size of a voxel along X-axis [m]
 */
#define VSX VSA
#endif

#ifndef VSY
/* *** Physical size of a voxel along Y-axis (up) [m]
 */
#define VSY VSA
#endif

#ifndef VSZ
/* *** Physical size of a voxel along Z-axis [m]
 */
#define VSZ VSA
#endif

#ifndef QRAD
/* *** Total radiation for a segment [MJ]
   The total radiation falling upon a segment. The actual intercepted radiation used
   for calculation of the physiological conditions is less than this and determined
   by the shading. Variation: 5.0-100.0, should be tested, no direct evidence,
   newly introduced parameter in LPFG-LIGNUM.
 */
#define QRAD 20.0
#endif

#ifndef SH
/* *** Shadow value for a segment [MJ]
   The shadow value formally is anti-radiation and, thus, has the units thereof. It
   is an absolute value, that is if SH = QRAD the self-shading, the shadow a segment
   exerts on itself, equals its maximum, i.e. the intercepted radiation is
   0. Variation: the same as for QRAD, but depends on its value, newly introduced
   parameter in LPFG-LIGNUM.
 */
#define SH 10.0
#endif

#ifndef SHLEN
/* *** Shadow length [m]
   The distance at which the shadow propagates. The shadow propagation model uses the
   pyramidal shadow propagation model from the exerting segment down through the
   voxel layers. The number of layers is approximated by SHLEN. The accuracy of the
   approximation also depends on the actual voxel size VSA. Variation: 0.3 - 1.0 is
   quite good choice, should be tested, newly introduced parameter in LPFG-LIGNUM.
 */
#define SHLEN 0.7
#endif

#ifndef NPOINCYL
/* *** Number of shadow sensors on a segment [int]
   Number of shadow sensors on a segment is used for more accurate assessment of the
   shadow conditions of a segment. Namely, the shadow exerted upon a segment is an
   average of the shadow values in each of the voxels that contain each of the
   sensors. Variation: 1 is enough.
*/
#define NPOINCYL 1
#endif

#ifndef LABELFLAG
/* *** Labels flag [int]
   Whether to draw the labels on the final plot: height, 1.0m-bar, diameter at
   base. 0 - NO, 1 (any non-zero) - YES
 */
#define LABELFLAG 1
#endif

#ifndef FOLFLAG
/* *** Foliage plotting flag [int]
   Whether to plot foliage. 0 - NO, 1 (any non-zero) - YES
 */
#define FOLFLAG 1
#endif
