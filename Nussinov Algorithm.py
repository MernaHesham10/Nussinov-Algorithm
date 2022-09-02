import numpy as npObject

def GetMatchedAminoAcid(AminoAcidInputToCheck):

    # Amino Acid Pair Dictionary Without Wobble Pairing
    AminoAcidPairDictionary = {
        'A':'U',
        'U':'A',
        'G':'C',
        'C':'G'
    }

    # Check If Input Amino Acid in AminoAcidPairDictionary
    if AminoAcidInputToCheck in AminoAcidPairDictionary:

        return AminoAcidPairDictionary[AminoAcidInputToCheck]

    return 0


def CreateNewSquareMatrix(RNAInputSequence):

    # Create Empty Matrix Called CreatedNewSquareMatrix Means Put nan In Each Cell in CreatedNewSquareMatrix
    CreatedNewSquareMatrix = npObject.empty([len(RNAInputSequence), len(RNAInputSequence)])
    CreatedNewSquareMatrix[:] = npObject.NAN

    # Put 0 in diagonal of CreatedNewSquareMatrix
    for RNAInputSequenceRowIndex in range(len(RNAInputSequence)):
        for RNAInputSequenceColIndex in range(len(RNAInputSequence)):
            if RNAInputSequenceColIndex == RNAInputSequenceRowIndex:
                CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] = 0

    # Put 0 in SubDiagonal of CreatedNewSquareMatrix
    for RNAInputSequenceRowIndex in range(1, len(RNAInputSequence)):
        for RNAInputSequenceColIndex in range(len(RNAInputSequence)):
            if RNAInputSequenceColIndex == RNAInputSequenceRowIndex - 1:
                CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] = 0

    return CreatedNewSquareMatrix


def CheckMatchedAminoAcid(AminoAcidInputToCheck):

    AminoAcidPairDictionary = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G"
    }

    if AminoAcidInputToCheck in AminoAcidPairDictionary.items():
        return True

    return False


def FillingSquareMatrixDiagonally(CreatedNewSquareMatrix, RNAInputSequence):

    #Create Bifurcation Dictionary to Store Max Value and Its Row, and Col
    bifurcationDictionary = {}
    for bifurcationDictionaryRowIndex in range(len(RNAInputSequence)):
        for bifurcationDictionaryColIndex in range(len(RNAInputSequence)):
            bifurcationDictionary[bifurcationDictionaryRowIndex, bifurcationDictionaryColIndex] = 0

    maxValuebOfBifurcation = 0
    MainLoopIndex = 0
    SquareMatrixRowIndex = 0
    SquareMatrixColIndex = 0
    bifurcationIndex = 0
    y = 1
    pairFound = 0

    #Get Max Value to Put in Matrix
    for MainLoopIndex in range(len(RNAInputSequence)):
        for SquareMatrixRowIndex in range(len(RNAInputSequence)):
            for SquareMatrixColIndex in range(1, len(RNAInputSequence)):

                if SquareMatrixColIndex == SquareMatrixRowIndex + y:

                    if GetMatchedAminoAcid(RNAInputSequence[SquareMatrixRowIndex]) == RNAInputSequence[SquareMatrixColIndex]:
                        pairFound = 1

                    elif ((SquareMatrixRowIndex < SquareMatrixColIndex) & (SquareMatrixRowIndex + 1 != SquareMatrixColIndex)):

                        for bifurcationIndex in range(SquareMatrixRowIndex + 1, SquareMatrixColIndex-1):

                            if maxValuebOfBifurcation < (CreatedNewSquareMatrix[SquareMatrixRowIndex][bifurcationIndex] + CreatedNewSquareMatrix[bifurcationIndex + 1][SquareMatrixColIndex]):

                                bifurcationDictionary[SquareMatrixRowIndex, SquareMatrixColIndex] = (CreatedNewSquareMatrix[SquareMatrixRowIndex][bifurcationIndex] + CreatedNewSquareMatrix[bifurcationIndex + 1][SquareMatrixColIndex])
                                maxValuebOfBifurcation = CreatedNewSquareMatrix[SquareMatrixRowIndex][bifurcationIndex] + CreatedNewSquareMatrix[bifurcationIndex + 1][SquareMatrixColIndex]

                    else:
                        pairFound = 0

                    CreatedNewSquareMatrix[SquareMatrixRowIndex][SquareMatrixColIndex] = max(CreatedNewSquareMatrix[SquareMatrixRowIndex][SquareMatrixColIndex - 1], CreatedNewSquareMatrix[SquareMatrixRowIndex + 1][SquareMatrixColIndex], CreatedNewSquareMatrix[SquareMatrixRowIndex + 1][SquareMatrixColIndex - 1] + pairFound, bifurcationDictionary[SquareMatrixRowIndex, SquareMatrixColIndex])
                    pairFound = 0
        y += 1

    print("CreatedNewSquareMatrix = \n", npObject.matrix(CreatedNewSquareMatrix))
    return CreatedNewSquareMatrix


def CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput, RNAInputSequenceRowIndex, RNAInputSequenceColIndex):
    bifurcationIndex = 0

    if RNAInputSequenceRowIndex < RNAInputSequenceColIndex:
        if CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] == CreatedNewSquareMatrix[RNAInputSequenceRowIndex + 1][RNAInputSequenceColIndex]: # 1st rule
            CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput, RNAInputSequenceRowIndex + 1, RNAInputSequenceColIndex)

        elif CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] == CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex - 1]: # 2nd rule
            CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput, RNAInputSequenceRowIndex, RNAInputSequenceColIndex - 1)


        elif CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] == CreatedNewSquareMatrix[RNAInputSequenceRowIndex + 1][RNAInputSequenceColIndex - 1] + CheckMatchedAminoAcid((RNAInputSequence[RNAInputSequenceRowIndex], RNAInputSequence[RNAInputSequenceColIndex])): # 3rd rule
            CalculatedTracebackOutput.append((RNAInputSequenceRowIndex, RNAInputSequenceColIndex))
            CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput, RNAInputSequenceRowIndex + 1, RNAInputSequenceColIndex - 1)

        else:
            for bifurcationIndex in range(RNAInputSequenceRowIndex + 1, RNAInputSequenceColIndex - 1):
                if CreatedNewSquareMatrix[RNAInputSequenceRowIndex][RNAInputSequenceColIndex] == CreatedNewSquareMatrix[RNAInputSequenceRowIndex,  bifurcationIndex] + CreatedNewSquareMatrix[ bifurcationIndex + 1][RNAInputSequenceColIndex]: # 4th rule
                    CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput, RNAInputSequenceRowIndex,  bifurcationIndex)
                    CalculateTracebackForMatrix(CreatedNewSquareMatrix, RNAInputSequence, CalculatedTracebackOutput,  bifurcationIndex + 1, RNAInputSequenceColIndex)
                    break
    return CalculatedTracebackOutput


def GetMatchedSymbol(RNAInputSequence, tracebackFunctionOutputList):
    dotCreatedList = ["." for i in range(len(RNAInputSequence))]

    for tracebackFunctionOutputListIndex in tracebackFunctionOutputList:
        dotCreatedList[min(tracebackFunctionOutputListIndex)] = "("
        dotCreatedList[max(tracebackFunctionOutputListIndex)] = ")"

    return "".join(dotCreatedList)


if __name__ == "__main__":

    # CGGACCCAGACUUUC
    # GGGAAAUCC
    RNAInputSequence = input("Please Enter The RNA Sequence: ")

    FillingSquareMatrixDiagonallyOut = FillingSquareMatrixDiagonally(CreateNewSquareMatrix(RNAInputSequence), RNAInputSequence)

    CalculatedTracebackOutput = []

    print(GetMatchedSymbol(RNAInputSequence, CalculateTracebackForMatrix(FillingSquareMatrixDiagonallyOut, RNAInputSequence, CalculatedTracebackOutput, 0, len(RNAInputSequence) - 1)))