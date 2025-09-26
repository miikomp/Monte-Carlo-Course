#define MT_IS_ELASTIC_SCATTER(mt)   ((mt) == 2)
#define MT_IS_FISSION(mt)           ((((mt) >= 18) && ((mt) <= 21)) || ((mt) == 38))
#define MT_IS_INELASTIC_SCATTER(mt) (((mt) >= 50) && ((mt) <= 91))
#define MT_IS_CAPTURE(mt)           (((mt) >= 102) && ((mt) <= 117))

#define MT_INELASTIC_CONTINUUM 91