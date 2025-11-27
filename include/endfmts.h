#define MT_IS_ELASTIC_SCATTER(mt)   ((mt) == 2)
#define MT_IS_FISSION(mt)           ((((mt) >= 18) && ((mt) <= 21)) || ((mt) == 38))
#define MT_IS_INELASTIC_SCATTER(mt) (((mt) >= 50) && ((mt) <= 91))
#define MT_IS_CAPTURE(mt)           (((mt) >= 102) && ((mt) <= 117))

#define MT_INELASTIC_CONTINUUM 91

/* Macroscopic detector response identifiers (negative values)
 * These group together related MTs instead of referring to a single MT.
 * -1: total
 * -2: total capture (MT_IS_CAPTURE)
 * -3: total elastic scatter (MT_IS_ELASTIC_SCATTER)
 * -4: total inelastic scatter (MT_IS_INELASTIC_SCATTER)
 * -5: total fission (MT_IS_FISSION)
 * Add new grouped macroscopic responses by defining an ID here and
 * extending the matcher table in computedetectorbin.c.
 */
#define DETECTOR_RESPONSE_MACRO_TOTAL        (-1)
#define DETECTOR_RESPONSE_MACRO_CAPTURE      (-2)
#define DETECTOR_RESPONSE_MACRO_ELASTIC      (-3)
#define DETECTOR_RESPONSE_MACRO_INELASTIC    (-4)
#define DETECTOR_RESPONSE_MACRO_FISSION      (-5)
