package column

import (
	"PainTheMaster/mybraly/mymath/round"
	"fmt"
	"math"
	"os"
)

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

//Init initializes a set of analyte+column
func (aic *AnalyteInCol) Init(
	iniConc float64,
	colLen float64,
	division int,
	flowRate float64,
	relRate float64,
	diffRate float64,
	decompRate float64,
	totalT float64,
	timeFrac float64,
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
