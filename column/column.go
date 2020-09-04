package column

import (
	"PainTheMaster/mybraly/mymath/round"
	"fmt"
	"math"
	"os"
)

//AnalyteInCol is a set of parametes--time-related configurations, analytical conditions and chemical properties-- and calcluation results.
type AnalyteInCol struct {
	//column: Primary
	column   []float64
	colLen   float64
	division int

	plateHight float64
	idxColHead int
	idxColTail int

	//Configuration: Primary
	flowRate float64
	relRate  float64

	//Analyte: Primary
	diffCoeff  float64
	decompRate float64

	//Analyte: configurational
	differential []float64
	idxDiffBegin int
	idxDiffEnd   int

	//Time
	totalT      float64
	totalCycl   int
	timeFrac    float64
	transCycl   int
	totalTransl int
}

//Init initializes a set of analyte+column. SI unit.
func (aic *AnalyteInCol) Init(
	iniConc float64, //M unit
	colLen float64, //m unit
	division int, //no dimension
	flowRate float64, //m/s
	relRate float64, //no dimension
	diffRate float64, //10^(-1) m^2 s^(-2)
	decompRate float64, //s^(-1)
	totalT float64, //s
	timeFrac float64, //s
) {
	aic.colLen = colLen
	aic.division = division
	aic.plateHight = aic.colLen / float64(aic.division)
	aic.idxColHead = division - 1
	aic.idxColTail = 0
	aic.differential = make([]float64, aic.division)

	aic.flowRate = flowRate
	aic.relRate = relRate

	aic.diffCoeff = diffRate
	aic.decompRate = decompRate

	aic.totalT = totalT
	aic.timeFrac = timeFrac

	aic.totalCycl = int(totalT / timeFrac)
	//	transCycl*dT*flowRate*relRate = plateHight
	aic.transCycl = int(round.Round(aic.plateHight / (timeFrac * flowRate * aic.relRate)))
	aic.totalTransl = aic.totalCycl / aic.transCycl
	aic.column = make([]float64, aic.totalTransl+aic.division)
	aic.column[aic.division-1] = iniConc
	aic.idxDiffBegin = aic.division - 2 //diffusable reasion
	aic.idxDiffEnd = aic.division - 1   //diffusable reasion
}

//Transl is translation.
func (aic *AnalyteInCol) Transl(cycl int) {
	if cycl%aic.transCycl == 0 {
		aic.idxColHead++
		aic.idxColTail++
		if aic.idxDiffBegin < aic.idxColTail {
			aic.idxDiffBegin = aic.idxColTail
		}
	}
}

//Diffuse is diffusion
func (aic *AnalyteInCol) Diffuse() {
	for i := aic.idxDiffBegin - aic.idxColTail; i <= aic.idxDiffEnd-1-aic.idxColTail; i++ {
		aic.differential[i] = (aic.column[i+1+aic.idxColTail] - aic.column[i+aic.idxColTail]) / aic.plateHight
	}
	for i := aic.idxDiffBegin - aic.idxColTail; i <= aic.idxDiffEnd-1-aic.idxColTail; i++ {
		flux := aic.differential[i] * aic.diffCoeff * aic.timeFrac
		aic.column[i+aic.idxColTail] += flux
		aic.column[i+1+aic.idxColTail] -= flux
	}
	if aic.idxDiffBegin >= aic.idxColTail+1 {
		aic.idxDiffBegin--
	}
	if aic.idxDiffEnd <= aic.idxColHead-1 {
		aic.idxDiffEnd++
	}
}

//React calculates reaction
func (aicThis *AnalyteInCol) React(aicThat *AnalyteInCol) {
	for i := 0; i <= aicThis.idxColHead-aicThis.idxColTail; i++ {
		decay := aicThis.column[aicThis.idxColTail+i] * (1.0 - math.Exp(-1.0*aicThis.decompRate*aicThis.timeFrac))
		aicThis.column[aicThis.idxColTail+i] -= decay
		aicThat.column[aicThat.idxColTail+i] += decay
	}
}

//Merge merges two chromatograms
func Merge(a, b AnalyteInCol) (merged AnalyteInCol) {
	var coarse, fine AnalyteInCol
	if a.timeFrac*float64(a.transCycl) >= b.timeFrac*float64(b.transCycl) {
		coarse = a
		fine = b
	} else {
		coarse = b
		fine = a
	}
	merged.division = coarse.division
	merged.timeFrac = coarse.timeFrac
	merged.transCycl = coarse.transCycl
	merged.totalTransl = coarse.totalTransl
	merged.column = append(coarse.column, []float64{}...)
	merged.idxColHead = coarse.idxColHead
	merged.idxColTail = coarse.idxColTail
	merged.idxDiffBegin = coarse.idxDiffBegin
	merged.idxDiffEnd = coarse.idxDiffEnd

	for i := 0; i <= merged.totalTransl-1; i++ {
		timeCoarse := float64(i*merged.transCycl) * merged.timeFrac
		idxFin := int(timeCoarse / (fine.timeFrac * float64(fine.transCycl)))
		left := float64(idxFin*fine.transCycl) * fine.timeFrac
		right := left + fine.timeFrac*float64(fine.transCycl)
		merged.column[i] += (fine.column[idxFin]*(right-timeCoarse) + fine.column[idxFin+1]*(timeCoarse-left)) / (fine.timeFrac * float64(fine.transCycl))
	}

	for i := 0; i <= merged.division-1; i++ {
		merged.column[merged.idxColTail+i] += coarse.column[coarse.idxColTail+i]
	}

	return
}

//Output outputs
func (aic AnalyteInCol) Output(fileName string) {
	file, err := os.Create(fileName)
	if err != nil {
		fmt.Println(err)
	}
	for i := 0; i <= aic.idxColTail-1; i++ {
		file.WriteString(fmt.Sprintf("%e,%e\n", aic.timeFrac*float64(aic.transCycl*i), aic.column[i]))
	}
	file.Close()
}

//RecomTimeDiv calcluatates recommended time division based on total T, colLen, division, flowRate, ralRate
func RecomTimeDiv(colLen float64, //m unit
	division int, //no dimension
	flowRate float64, //m/s
	relRate float64, //no dimension
	totalT float64, //s
) (optTimeDiv int) {
	const optTranslCycl = 100

	cellLength := colLen / float64(division)
	optTimeFrac := cellLength / (optTranslCycl * flowRate * relRate)
	optTimeDiv = int(totalT / optTimeFrac)
	return
}
