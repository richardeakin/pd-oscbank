// ------------------------ oscbank~ 0.2 -----------------------------
// oscillator bank using 3 seperate float inlets with interpolation
// - Please see the included README and LICENSE for more info
// - reakinator@gmail.com
// -------------------------------------------------------------------

#include "m_pd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef NT
#pragma warning( disable : 4244 )
#define inline
#endif

#define WAVETABLESIZE 65536 // 2^16
#define DEFAULT_NPARTIALS 100
#define DEFAULT_interp_incr 0.0045 // per sample, this is 20 ms @ 44k sr
#define NEAR_ZERO 0.0000001

typedef struct _partial
{
	int   index;
	float fCurr;
	float freq;
	float fIncr;
	float aCurr;
	float amp;
	float aIncr;
	float phase;
	unsigned long  nInterp;
} t_partial;

typedef struct _oscbank
{
	t_object x_obj;
	t_outlet *list_outlet;
	float    *wavetable;
	int      wavetablesize;
	int      got_a_table;
	t_partial *pBank;
	float    infreq;
	float    inamp;
	float    sampleRate;
	float    sampleperiod;
	float    interp_incr;
	long     interpSamples;
	int      nPartials;
} t_oscbank;

static t_class *oscbank_class;
static void *oscbank_new(void);
static void oscbank_free(t_oscbank *x);
static void oscbank_interpMs(t_oscbank *x, t_floatarg n);
static void oscbank_nPartials(t_oscbank *x, t_floatarg n);
static void oscbank_index(t_oscbank *x, t_floatarg in);
static void oscbank_table(t_oscbank *x, t_symbol *tablename);
static void oscbank_print(t_oscbank *x);
static void oscbank_outlist(t_oscbank *x);
static void oscbank_dsp(t_oscbank *x, t_signal **sp);
static void oscbank_reset(t_oscbank *x);

// ------------------- Setup / Teardown ------------------------------

void oscbank_tilde_setup(void)
{
	oscbank_class = class_new(gensym("oscbank~"), (t_newmethod)oscbank_new, (t_method)oscbank_free,
			sizeof(t_oscbank), 0, A_DEFFLOAT, 0);
	class_addfloat(oscbank_class, oscbank_index);
	class_addmethod(oscbank_class, (t_method)oscbank_table, gensym("table"), A_SYMBOL);
	class_addmethod(oscbank_class, (t_method)oscbank_interpMs, gensym("interp"), A_FLOAT, 0);
	class_addmethod(oscbank_class, (t_method)oscbank_dsp, gensym("dsp"), (t_atomtype)0);
	class_addmethod(oscbank_class, (t_method)oscbank_print, gensym("print"), 0);
	class_addmethod(oscbank_class, (t_method)oscbank_outlist, gensym("sendout"), 0);
	class_addmethod(oscbank_class, (t_method)oscbank_reset, gensym("reset"), 0);
	class_addmethod(oscbank_class, (t_method)oscbank_nPartials, gensym("partials"), A_FLOAT, 0);
}

static void *oscbank_new(void)
{
	t_oscbank *x = (t_oscbank *)pd_new(oscbank_class);

	float twopi;
	int i;

	outlet_new(&x->x_obj, gensym("signal"));
	x->list_outlet = outlet_new(&x->x_obj, gensym("list"));

	floatinlet_new(&x->x_obj, &x->infreq);
	floatinlet_new(&x->x_obj, &x->inamp);
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("interp"));

	// hardcoded value prevents devide by zero in oscbank_index(), but will be updated when ddsp is switched on
	x->sampleRate = 44100;
	x->sampleperiod = 1.0f / x->sampleRate;
	oscbank_interpMs( x, 20.0f);

	x->got_a_table = 0;
	x->nPartials = DEFAULT_NPARTIALS;
	x->pBank = (t_partial *)getbytes( x->nPartials * sizeof(t_partial));
	memset(x->pBank, 0, x->nPartials * sizeof(t_partial));

	twopi = 8.0f * atan(1.0f);
	x->wavetablesize = WAVETABLESIZE;
	float *sinewave;
	sinewave = (t_float *)malloc(x->wavetablesize * sizeof(t_float));
	for(i = 0; i < x->wavetablesize; i++)
		sinewave[i] = sin(twopi * (float)i/ x->wavetablesize);

	x->wavetable = &sinewave[0];

	return (x);
}

static void oscbank_free(t_oscbank *x)
{
	free(x->pBank);
	if(!x->got_a_table) {
		free(x->wavetable);
	}
}

// ------------------- External Message Routines --------------------

/*
 * Interpolation Time:
 * milleseconds to interpolate over; so samples = (n*SR)/1000
 * divide only when converting the interp time to samples(here),since it 
 * is only used as a denominator to find the increment proportion:
 * SP= 1/SR, 1/(n*SR/1000) = (1000*SP)/n
 */
static void oscbank_interpMs(t_oscbank *x, t_floatarg n)
{
	if(n > 0) {
		x->interp_incr = (1000 * x->sampleperiod) / n ;	
	} else {
		x->interp_incr = x->sampleperiod; 	
	}
	x->interpSamples = (unsigned long)((n *.001) * x->sampleRate);
}

static void oscbank_nPartials(t_oscbank *x, t_floatarg n)
{
	x->pBank = (t_partial *)resizebytes( x->pBank, x->nPartials * sizeof(t_partial), n * sizeof(t_partial));
	x->nPartials = n;
	post("max partials: %d", x->nPartials);
}

static void oscbank_index(t_oscbank *x, t_floatarg in)
{
	int i, iindex;
	iindex = (int)in;
	t_partial *bank = x->pBank;
	int empty_index = -1;
	int quietest_index = 0;

	if( iindex < 0)	{
		error("oscbank~ needs a positive index.");
		return;
	}

	//check if it is continuing partial
	//recaluclate increment slope from current interpolated positions and update goal
	for(i = 0; i < x->nPartials; i++) {
		if(bank[i].index == iindex) {
			if(bank[i].aCurr == 0) bank[i].aCurr = NEAR_ZERO;
			bank[i].fIncr = (x->infreq - bank[i].fCurr) * x->interp_incr;
			bank[i].aIncr = (x->inamp - bank[i].aCurr) * x->interp_incr;
			bank[i].freq = x->infreq;
			bank[i].amp = x->inamp;
			bank[i].nInterp = x->interpSamples;
			return;
		}
		// store indeces either for an empty partial or quietest
		if(empty_index < 0 && bank[i].aCurr < NEAR_ZERO) {
			empty_index = i;
		}
		if(bank[i].amp < bank[quietest_index].amp) {
			quietest_index = i;
		}
	}

	if(empty_index >= 0) {
		//ramp amp from zero
		bank[empty_index].index = iindex;
		bank[empty_index].fCurr = x->infreq;
		bank[empty_index].fIncr = 0;
		bank[empty_index].freq = x->infreq;
		bank[empty_index].amp = x->inamp;
		bank[empty_index].nInterp = x->interpSamples;
		bank[empty_index].aCurr = NEAR_ZERO;  
		bank[empty_index].aIncr = x->inamp * x->interp_incr;
		return;
	}

	//oscbank is full, steal quietest partial and ramp amp from zero
	bank[quietest_index].index = iindex;
	bank[quietest_index].fCurr = x->infreq;
	bank[quietest_index].fIncr = 0;
	bank[quietest_index].freq = x->infreq;
	bank[quietest_index].amp = x->inamp;
	bank[quietest_index].nInterp = x->interpSamples;
	bank[quietest_index].aCurr = NEAR_ZERO;
	bank[quietest_index].aIncr = x->inamp * x->interp_incr; 
}

static void oscbank_table(t_oscbank *x, t_symbol *tablename)
{
	if(!x->got_a_table) {
		free(x->wavetable);
		x->got_a_table = 0;
	}

	t_garray *a;
	if (!(a = (t_garray *)pd_findbyclass(tablename, garray_class))) {
		pd_error(x, "%s: no such array", tablename->s_name);
	} else if (!garray_getfloatarray(a, &x->wavetablesize, &x->wavetable)) {
		pd_error(x, "%s: bad template for tabread", tablename->s_name);
	} else {
		post("wavetablesize: %d", x->wavetablesize );
		x->got_a_table = 1;
	}
}

static void oscbank_print(t_oscbank *x)
{
	t_partial *bank = x->pBank;
	post("#:  Index,  Freq,  Amp");

	int i;
	for(i = 0; i < x->nPartials; i++) {
		if(bank[i].aCurr) {
			post("%d: %d, %.2f, %.2f", i, bank[i].index, bank[i].freq,  bank[i].amp );
		}
	}
}

static void oscbank_outlist(t_oscbank *x)
{
	t_partial *bank = x->pBank;
	t_atom *outv;
	int audiblePartials = 0;
	outv = (t_atom *)malloc(x->nPartials * 3 * sizeof(t_atom));
	if(!outv) {
		pd_error(x, "could not allocate memory for t_atom list");
		return;
	}

	int i, offset;
	for(i = 0; i < x->nPartials; i++) {
		if(bank[i].aCurr) {
			offset = audiblePartials++ * 3;
			SETFLOAT(outv + offset, (float)bank[i].index);
			SETFLOAT(outv + offset + 1, bank[i].freq);
			SETFLOAT(outv + offset + 2, bank[i].amp);
		}
	}
	outlet_list(x->list_outlet, &s_list, audiblePartials * 3, outv); // only send out atoms that are filled
	free(outv);
}

static void oscbank_reset(t_oscbank *x)
{
	memset(x->pBank, 0, x->nPartials * sizeof(t_partial));
}

//------------------------- DSP routines ---------------------------------------

static t_int *oscbank_perform(t_int *w)
{
	t_oscbank *x = (t_oscbank *)(w[1]);
	t_float *out = (t_float *)(w[2]);
	t_int n = (t_int)(w[3]);
	t_int i, sample;
	t_float phaseincrement;
	t_int	lookup;
	t_partial *bank = x->pBank;

	memset(out, 0, n *sizeof( t_float ));

	for(i=0; i < x->nPartials; i++) {
		if(bank[i].aCurr != 0) {
			for(sample = 0; sample < n; sample++) {
				if(bank[i].nInterp > 0) {
					bank[i].fCurr += bank[i].fIncr;
					bank[i].aCurr += bank[i].aIncr;
					--bank[i].nInterp;
				} else {
					bank[i].fCurr = bank[i].freq;
					bank[i].aCurr = bank[i].amp;
				}

				// get the phase increment freq = cyc/sec,
				// sr = samp/sec, phaseinc = cyc/samp = freq/sr = freq * sampleperiod
				phaseincrement = bank[i].fCurr * x->sampleperiod;
				bank[i].phase += phaseincrement;
				while(bank[i].phase >= 1.0f) {
					bank[i].phase -= 1.0f;	
				}
				while(bank[i].phase < 0.0f) {
					bank[i].phase += 1.0f;
				}
				lookup = (int)(x->wavetablesize * bank[i].phase); 
				*(out+sample) += *(x->wavetable + lookup) * bank[i].aCurr; 
			}
		}
	}
	return (w+4);
}

static void oscbank_dsp(t_oscbank *x, t_signal **sp)
{
	x->sampleRate =  sp[0]->s_sr;
	x->sampleperiod = 1 / x->sampleRate;
	dsp_add(oscbank_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}
