NULL


AAlist <- c("A", "C", "D", "E", "F", "G", "H",
    "I", "K", "L", "M", "N", "P", "Q", "R", "S",
    "T", "V", "W", "Y")
# ----------------------------------------------------------------------------
# # MHC-II contacting residue groups for DRB
# molecules
p1 <- list(c(111, 114, 115, 118), c(111, 114,
    115, 118), c(111, 114, 115, 118), c(111, 114,
    115, 118))
G1 <- list(c("N", "A", "V", "F"), c("Y", "G",
    "E", "T"), c("N", "V", "G", "F"), c("N", "V",
    "V", "F"))
p2 <- list(c(106, 107, 110, 111), c(106, 107,
    110, 111), c(106, 107, 110, 111), c(106, 107,
    110, 111), c(106, 107, 110, 111))
G2 <- list(c("T", "Y", "H", "N"), c("N", "Y",
    "H", "N"), c("T", "V", "H", "N"), c("Y", "C",
    "N", "Y"), c("T", "Y", "Y", "N"))
p3 <- list(c(103, 107), c(103, 107), c(103, 107),
    c(103, 107), c(103, 107), c(103, 107), c(103,
        107), c(103, 107))
G3 <- list(c("A", "Y"), c("E", "V"), c("E", "Y"),
    c("L", "Y"), c("Q", "V"), c("Q", "Y"), c("R",
        "Y"), c("V", "C"))
p4 <- list(c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107), c(40, 42, 43, 55, 57, 99, 100, 103,
    107))

G4 <- list(c("L", "F", "E", "L", "E", "Q", "R",
    "A", "Y"), c("S", "S", "E", "Y", "D", "Q",
    "K", "R", "Y"), c("S", "S", "E", "F", "E",
    "Q", "K", "R", "Y"), c("V", "H", "E", "F",
    "D", "Q", "K", "A", "Y"), c("V", "H", "E",
    "F", "D", "Q", "R", "A", "Y"), c("G", "Y",
    "K", "F", "E", "D", "R", "Q", "V"), c("S",
    "G", "E", "F", "D", "D", "R", "L", "Y"), c("S",
    "G", "E", "F", "D", "R", "R", "V", "C"), c("D",
    "F", "E", "Y", "H", "R", "R", "E", "V"), c("S",
    "S", "E", "F", "D", "D", "R", "A", "Y"), c("S",
    "G", "E", "L", "E", "D", "R", "A", "Y"), c("S",
    "S", "E", "F", "D", "D", "E", "A", "Y"), c("S",
    "S", "E", "F", "E", "Q", "R", "A", "Y"), c("S",
    "G", "E", "F", "D", "R", "R", "E", "Y"), c("S",
    "S", "E", "F", "E", "D", "R", "L", "Y"), c("P",
    "R", "E", "F", "D", "Q", "A", "A", "Y"), c("R",
    "S", "E", "Y", "D", "Q", "K", "R", "Y"), c("L",
    "S", "E", "F", "E", "Q", "K", "Q", "Y"), c("A",
    "C", "E", "N", "I", "R", "R", "E", "Y"), c("D",
    "Y", "E", "F", "H", "D", "R", "A", "Y"))
p5 <- list(c(40, 42, 57, 59, 99, 100, 103), c(40,
    42, 57, 59, 99, 100, 103), c(40, 42, 57, 59,
    99, 100, 103), c(40, 42, 57, 59, 99, 100,
    103), c(40, 42, 57, 59, 99, 100, 103), c(40,
    42, 57, 59, 99, 100, 103), c(40, 42, 57, 59,
    99, 100, 103), c(40, 42, 57, 59, 99, 100,
    103), c(40, 42, 57, 59, 99, 100, 103), c(40,
    42, 57, 59, 99, 100, 103), c(40, 42, 57, 59,
    99, 100, 103), c(40, 42, 57, 59, 99, 100,
    103), c(40, 42, 57, 59, 99, 100, 103), c(40,
    42, 57, 59, 99, 100, 103), c(40, 42, 57, 59,
    99, 100, 103), c(40, 42, 57, 59, 99, 100,
    103), c(40, 42, 57, 59, 99, 100, 103), c(40,
    42, 57, 59, 99, 100, 103), c(40, 42, 57, 59,
    99, 100, 103), c(40, 42, 57, 59, 99, 100,
    103))
G5 <- list(c("L", "F", "E", "C", "Q", "R", "A"),
    c("S", "S", "D", "Y", "Q", "K", "R"), c("S",
        "S", "E", "Y", "Q", "K", "R"), c("V",
        "H", "D", "Y", "Q", "K", "A"), c("V",
        "H", "D", "Y", "Q", "R", "A"), c("G",
        "Y", "E", "L", "D", "R", "Q"), c("S",
        "G", "D", "Y", "D", "R", "L"), c("S",
        "G", "D", "Y", "R", "R", "V"), c("D",
        "F", "H", "G", "R", "R", "E"), c("S",
        "S", "D", "Y", "D", "R", "A"), c("S",
        "G", "E", "H", "D", "R", "A"), c("S",
        "S", "D", "Y", "D", "E", "A"), c("S",
        "S", "E", "Y", "Q", "R", "A"), c("S",
        "G", "D", "Y", "R", "R", "E"), c("S",
        "S", "E", "Y", "D", "R", "L"), c("P",
        "R", "D", "Y", "Q", "A", "A"), c("R",
        "S", "D", "Y", "Q", "K", "R"), c("L",
        "S", "E", "Y", "Q", "K", "Q"), c("A",
        "C", "I", "Y", "R", "R", "E"), c("D",
        "Y", "H", "D", "D", "R", "A"))
p6 <- list(c(38, 40, 42, 57, 59, 99, 100, 13),
    c(38, 40, 42, 57, 59, 99, 100, 13), c(38,
        40, 42, 57, 59, 99, 100, 13), c(38, 40,
        42, 57, 59, 99, 100, 13), c(38, 40, 42,
        57, 59, 99, 100, 13), c(38, 40, 42, 57,
        59, 99, 100, 13), c(38, 40, 42, 57, 59,
        99, 100, 13), c(38, 40, 42, 57, 59, 99,
        100, 13), c(38, 40, 42, 57, 59, 99, 100,
        13), c(38, 40, 42, 57, 59, 99, 100, 13),
    c(38, 40, 42, 57, 59, 99, 100, 13), c(38,
        40, 42, 57, 59, 99, 100, 13), c(38, 40,
        42, 57, 59, 99, 100, 13), c(38, 40, 42,
        57, 59, 99, 100, 13), c(38, 40, 42, 57,
        59, 99, 100, 13), c(38, 40, 42, 57, 59,
        99, 100, 13), c(38, 40, 42, 57, 59, 99,
        100, 13), c(38, 40, 42, 57, 59, 99, 100,
        13), c(38, 40, 42, 57, 59, 99, 100, 13),
    c(38, 40, 42, 57, 59, 99, 100, 13))
G6 <- list(c("W", "L", "F", "E", "C", "Q", "R",
    "A"), c("E", "S", "S", "D", "Y", "Q", "K",
    "R"), c("E", "S", "S", "E", "Y", "Q", "K",
    "R"), c("E", "V", "H", "D", "Y", "Q", "K",
    "A"), c("E", "V", "H", "D", "Y", "Q", "R",
    "A"), c("W", "G", "Y", "E", "L", "D", "R",
    "Q"), c("E", "S", "G", "D", "Y", "D", "R",
    "L"), c("E", "S", "G", "D", "Y", "R", "R",
    "V"), c("K", "D", "F", "H", "G", "R", "R",
    "E"), c("E", "S", "S", "D", "Y", "D", "R",
    "A"), c("E", "S", "G", "E", "H", "D", "R",
    "A"), c("E", "S", "S", "D", "Y", "D", "E",
    "A"), c("E", "S", "S", "E", "Y", "Q", "R",
    "A"), c("E", "S", "G", "D", "Y", "R", "R",
    "E"), c("E", "S", "S", "E", "Y", "D", "R",
    "L"), c("W", "P", "R", "D", "Y", "Q", "A",
    "A"), c("E", "R", "S", "D", "Y", "Q", "K",
    "R"), c("E", "L", "S", "E", "Y", "Q", "K",
    "Q"), c("E", "A", "C", "I", "Y", "R", "R",
    "E"), c("Q", "D", "Y", "H", "D", "D", "R",
    "A"))

p7 <- list(c(40, 57, 59, 76, 90, 96, 99, 100),
    c(40, 57, 59, 76, 90, 96, 99, 100), c(40,
        57, 59, 76, 90, 96, 99, 100), c(40, 57,
        59, 76, 90, 96, 99, 100), c(40, 57, 59,
        76, 90, 96, 99, 100), c(40, 57, 59, 76,
        90, 96, 99, 100), c(40, 57, 59, 76, 90,
        96, 99, 100), c(40, 57, 59, 76, 90, 96,
        99, 100), c(40, 57, 59, 76, 90, 96, 99,
        100), c(40, 57, 59, 76, 90, 96, 99, 100),
    c(40, 57, 59, 76, 90, 96, 99, 100), c(40,
        57, 59, 76, 90, 96, 99, 100), c(40, 57,
        59, 76, 90, 96, 99, 100), c(40, 57, 59,
        76, 90, 96, 99, 100), c(40, 57, 59, 76,
        90, 96, 99, 100), c(40, 57, 59, 76, 90,
        96, 99, 100), c(40, 57, 59, 76, 90, 96,
        99, 100), c(40, 57, 59, 76, 90, 96, 99,
        100), c(40, 57, 59, 76, 90, 96, 99, 100),
    c(40, 57, 59, 76, 90, 96, 99, 100), c(40,
        57, 59, 76, 90, 96, 99, 100))
G7 <- list(c("L", "E", "C", "Y", "W", "L", "Q",
    "R"), c("S", "D", "Y", "F", "W", "L", "Q",
    "K"), c("S", "E", "Y", "Y", "W", "L", "Q",
    "K"), c("V", "D", "Y", "Y", "W", "L", "Q",
    "K"), c("V", "D", "Y", "Y", "W", "L", "Q",
    "R"), c("G", "E", "L", "Y", "W", "I", "D",
    "R"), c("S", "D", "Y", "Y", "W", "F", "D",
    "R"), c("S", "D", "Y", "Y", "W", "L", "D",
    "R"), c("S", "D", "Y", "Y", "W", "L", "R",
    "R"), c("D", "H", "G", "Y", "W", "F", "R",
    "R"), c("S", "D", "Y", "F", "W", "F", "D",
    "R"), c("S", "E", "H", "F", "W", "I", "D",
    "R"), c("S", "E", "H", "F", "W", "F", "D",
    "R"), c("S", "D", "Y", "F", "W", "I", "D",
    "E"), c("S", "E", "Y", "Y", "W", "L", "Q",
    "R"), c("S", "E", "Y", "Y", "W", "L", "D",
    "R"), c("P", "D", "Y", "F", "W", "I", "Q",
    "A"), c("R", "D", "Y", "Y", "W", "L", "Q",
    "K"), c("L", "E", "Y", "Y", "W", "L", "Q",
    "K"), c("A", "I", "Y", "Y", "W", "L", "R",
    "R"), c("D", "H", "D", "Y", "W", "F", "D",
    "R"))
p8 <- list(c(89, 90), c(89, 90), c(89, 90))
G8 <- list(c("H", "W"), c("S", "W"), c("Y", "W"))
p9 <- list(c(38, 59, 66, 67, 86, 89, 90), c(38,
    59, 66, 67, 86, 89, 90), c(38, 59, 66, 67,
    86, 89, 90), c(38, 59, 66, 67, 86, 89, 90),
    c(38, 59, 66, 67, 86, 89, 90), c(38, 59, 66,
        67, 86, 89, 90), c(38, 59, 66, 67, 86,
        89, 90), c(38, 59, 66, 67, 86, 89, 90),
    c(38, 59, 66, 67, 86, 89, 90), c(38, 59, 66,
        67, 86, 89, 90), c(38, 59, 66, 67, 86,
        89, 90), c(38, 59, 66, 67, 86, 89, 90),
    c(38, 59, 66, 67, 86, 89, 90), c(38, 59, 66,
        67, 86, 89, 90))
G9 <- list(c("W", "C", "S", "V", "D", "Y", "W"),
    c("E", "Y", "N", "V", "D", "Y", "W"), c("E",
        "Y", "Y", "V", "D", "Y", "W"), c("E",
        "Y", "Y", "V", "S", "Y", "W"), c("W",
        "L", "F", "V", "V", "S", "W"), c("E",
        "Y", "Y", "V", "I", "Y", "W"), c("K",
        "G", "N", "V", "V", "S", "W"), c("E",
        "H", "L", "L", "V", "S", "W"), c("E",
        "Y", "F", "V", "A", "H", "W"), c("W",
        "Y", "S", "V", "D", "Y", "W"), c("E",
        "Y", "F", "L", "V", "S", "W"), c("E",
        "Y", "F", "V", "V", "S", "W"), c("E",
        "Y", "Y", "A", "D", "Y", "W"), c("Q",
        "D", "D", "L", "D", "Y", "W"))
P <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9)
G <- list(G1, G2, G3, G4, G5, G6, G7, G8, G9)

# ----------------------------------------------------------------------------
# # MHC-II Polymorphic residue group for DP
# molecules
dp1 <- list(c(116, 119), c(116, 119))
dG1 <- list(c("V", "Q"), c("M", "Q"))
dp2 <- list(c(108))
dG2 <- list(c("H"))
dp3 <- list(c(108))
dG3 <- list(c("H"))
dp4 <- list(c(101, 108))
dG4 <- list(c("V", "H"))
dp5 <- list(c(101))
dG5 <- list(c("V"))
dp6 <- list(c(38, 100), c(38, 100))
dG6 <- list(c("Y", "A"), c("F", "A"))
dp7 <- list(c(100))
dG7 <- list(c("A"))
dp8 <- list()
dG8 <- list()  #off
dp9 <- list(c(38, 66, 67, 86), c(38, 66, 67, 86))
dG9 <- list(c("Y", "R", "F", "E"), c("F", "R",
    "F", "E"))
dP <- list(dp1, dp2, dp3, dp4, dp5, dp6, dp7,
    dp8, dp9)
dG <- list(dG1, dG2, dG3, dG4, dG5, dG6, dG7,
    dG8, dG9)
# ----------------------------------------------------------------------------
# # MHC-II polymorphic residue groups for DQ
# molecules
qp1 <- list(c(121, 124), c(121, 124))
qG1 <- list(c("T", "Q"), c("G", "Q"))
qp2 <- list(c(113))
qG2 <- list(c("H"))
qp3 <- list(c(113))
qG3 <- list(c("H"))
qp4 <- list(c(106, 113), c(106, 113), c(106, 113))
qG4 <- list(c("A", "H"), c("E", "H"), c("S", "H"))
qp5 <- list(c(106), c(106), c(106))
qG5 <- list(c("A"), c("E"), c("S"))
qp6 <- list(c(44, 106), c(44, 106), c(44, 106))
qG6 <- list(c("K", "A"), c("K", "E"), c("K", "S"))
qp7 <- list(c(106), c(106), c(106))
qG7 <- list(c("A"), c("E"), c("S"))
qp8 <- list()
qG8 <- list()  #off
qp9 <- list(c(44, 72, 73, 92))
qG9 <- list(c("K", "F", "D", "Y"))
qP <- list(qp1, qp2, qp3, qp4, qp5, qp6, qp7,
    qp8, qp9)
qG <- list(qG1, qG2, qG3, qG4, qG5, qG6, qG7,
    qG8, qG9)

#' Parameters (Hamiltonians) vector
#'
#' This function setups a vector containing Hamiltonians
#' that are corresponding to the given locus.
#' @param  locus One of the three loci: DQ, DP, or DRB.
#' @return A dataframe of two columns, the first column contains the names of
#'  the Hamiltonians and the second column contains their values. By default
#'  these values are set to zero.
#' @examples This is how to call this function
#'   parametersVector('DQ')
#'
#' @export
parametersVector <- function(locus) {
    locus <- toupper(locus)
    if (!(locus %in% c("DRB", "DP", "DQ"))) {
        stop(" Please choose the correct locus: DRB, DP, or DQ")
    }
    l <- c("b0")
    if (locus == "DRB") {
        for (i in AAlist) {
            l <- c(l, paste("H(", i, ")", sep = ""))
        }
        for (i in AAlist) {
            for (j in 1:length(P)) {
                for (k in 1:length(P[[j]])) {
                  l <- c(l, paste("H(", i, ",",
                    "g", j, "_", k, ")", sep = ""))
                }
            }
        }
    } else if (locus == "DP") {
        for (i in AAlist) {
            l <- c(l, paste("H(", i, ")", sep = ""))
        }
        for (i in AAlist) {
            for (j in 1:length(dP)) {
                for (k in 1:length(dP[[j]])) {
                  l <- c(l, paste("H(", i, ",",
                    "g", j, "_", k, ")", sep = ""))
                }
            }
        }

    } else {
        for (i in AAlist) {
            l <- c(l, paste("H(", i, ")", sep = ""))
        }
        for (i in AAlist) {
            for (j in 1:length(qP)) {
                for (k in 1:length(qP[[j]])) {
                  l <- c(l, paste("H(", i, ",",
                    "g", j, "_", k, ")", sep = ""))
                }
            }
        }
    }
    p <- rep(0, length(l))
    names(p) <- l
    return(as.data.frame(p))
}
# ----------------------------------------------------------------------------

#' Getting all possible frames
#'
#' This fucntion extracts all the possible frames  from the
#'  giveing peptide (frames of 9 mers).
#' @param peptide The amino acids sequence of  peptide of any length >= 9 mer.
#' @return  A vector containing all the possible frames on the given peptide.
#' @examples This is how to call this function
#'     sliding('HKLMNARTIPPL')
#'
#' @export
sliding <- function(peptide) {
    n <- nchar(peptide) - 8
    slides <- c()
    for (i in 1:n) {
        slides <- c(slides, substr(peptide, i,
            i + 8))
    }
    return(slides)
}


#' Second order Hamiltonians
#'
#' This function  processes the second-order Hamiltonians between the peptide
#' and mhc-ii molecule.
#'
#' @param peptide The peptide sequence
#' @param MHC The beta sequence of the targeted mhc-ii molecule
#' @param locus The intended mhc-ii locus, and it must be on of the
#' following loci: DQ, DP, and  DRB.
#' @return A vector conatining the  pairwise interacting Hamiltonians
#' @examples This is how to use this function
#'  secondOrderHamiltonians(pep, mhc, drb)
#'
#' @export
secondOrderHamiltonians <- function(peptide, MHC,
    locus) {
    peptide <- unlist(strsplit(toupper(peptide),
        split = ""))
    MHC <- unlist(strsplit(toupper(MHC), split = ""))
    L <- c()
    if (locus == "DRB") {
        for (i in 1:9) {
            for (j in 1:length(P[[i]])) {
                if (identical(MHC[P[[i]][[j]]],
                  G[[i]][[j]]) == T) {
                  L <- c(L, paste("H(", peptide[i],
                    ",", "g", i, "_", j, ")",
                    sep = ""))
                }
            }
        }
    } else if (locus == "DP") {
        for (i in 1:9) {
            if (i != 8) {
                for (j in 1:length(dP[[i]])) {
                  if (identical(MHC[dP[[i]][[j]]],
                    dG[[i]][[j]]) == T) {
                    L <- c(L, paste("H(", peptide[i],
                      ",", "g", i, "_", j, ")",
                      sep = ""))
                  }
                }
            }
        }
    } else {
        for (i in 1:9) {
            if (i != 8) {
                for (j in 1:length(qP[[i]])) {
                  if (identical(MHC[qP[[i]][[j]]],
                    qG[[i]][[j]]) == T) {
                    L <- c(L, paste("H(", peptide[i],
                      ",", "g", i, "_", j, ")",
                      sep = ""))
                  }
                }
            }
        }
    }
    return(L)
}


#' Binding vector
#'
#' This fucntion processes the  binding vector between the peptide and mhc-ii
#' molecule from all possible registers on the giving peptide.
#'
#' @param peptide The peptide in terms of amino acids sequence
#' @param MHC The beta sequence of the mhc-ii molecule
#' @param locue The type of mhc-ii molecule
#' @return  A dataframe contains both the first- and second-oreder
#'  Hamiltonians
#' @examples This is how to call this function
#' binding Vector('KLAMDVNNKHFRT', 'DRB1_0101', 'DRB')
#'
#' @export
bindingVector <- function(peptide, MHC, locus) {
    Pep <- unlist(strsplit(toupper(peptide), split = ""))
    mhC <- unlist(strsplit(toupper(MHC), split = ""))
    q <- c()
    for (j in 1:length(Pep)) {
        q <- c(q, paste("H(", Pep[j], ")", sep = ""))
    }

    slides <- sliding(peptide)
    for (s in 1:length(slides)) {
        q <- c(q, secondOrderHamiltonians(slides[s],
            MHC, locus))
    }
    p <- parametersVector(locus)
    for (r in rownames(p)) {
        p[r, ] <- length(which(q == r))
    }
    p["b0", 1] <- 1
    return(p)
}


#' Getting the actual binding vector
#'
#' This function creates the actual binding vector for peptide:MHC-II
#' interaction.
#' @return  A dataframe of two columns, the first column contains the names of
#'  the Hamiltonians and the second column contains their obtained values.
#' @param peptide The peptide sequences
#' @param mhcName The the name of the mhc-ii molecule
#' @param locue The type of mhc-ii molecule
#' @examples
#'  getInputVector('KHVNWQERTPLK', 'DRB1_0404_seq', 'DRB')
#'
#' @export
getInputVector <- function(peptide, mhcName, locus) {
    # dg <- system.file('extdata',
    # 'MHCII_Molecules_Sequences.Rda', package =
    # getPackageName())
    dg <- readRDS("./inst/extdata/MHCII_Molecules_Sequences.Rda")
    if (!(mhcName %in% dg[, 1])) {
        stop("Currently, the package does not support the binding specifities for this molecule.
            Please use the supportedMolecules to see the covered molecules")
    }
    v <- bindingVector(peptide, dg[which(dg$Pre_name ==
        mhcName), 3], locus)
    return(v)
}
# ----------------------------------------------------------------------------
#' Supported MHC class II molecules
#'
#' This function shows  MHC-II molecules  covered within this package.
#' @return  A dataframe containing two columns, the first column shows MHC-II
#' molecule name and the second column shows its alias neme. Please use
#' the names that appear in first column
#' @examples
#' ## What are the MHC-II molecules covered by this package?
#' supportedMolecules()
#'
#' @export
supportedMolecules <- function() {
    # dg <- system.file('extdata',
    # 'MHCII_Molecules_Sequences.Rda', package =
    # getPackageName())
    dg <- readRDS("./inst/extdata/MHCII_Molecules_Sequences.Rda")
    mhcs <- data.frame(dg[1], dg[2])
    colnames(mhcs) <- c("Supported Molecules",
        "Alias Name")
    return(mhcs)
}

#' Predicts peptide:MHC-II binding interaction.
#'
#' This function predicts the outcome of interaction between peptide
#' and mhc-ii
#' molecule by calculating the energy of the binding from all possible
#' configurations between the peptide and the mhc-ii ligand. It derives such
#' configurations by using function getInputVector
#' @param pep The amino acids sequence of the peptide of length =>9 nano mer.
#' @param mhc_name The name of mhc-ii molecule.
#' @param locus One of the three standard loci: DQ, DP, and DRB
#' @param opt This parameter determines the prediction type (intra-allele
#' or trans-alleles). 1 stands for intra-alleleic and 2 for trans-allelic.
#' The default value of this is 1.
#' @examples
#' predictBinding('AIGIITLYLGAVVQA', 'DRB1_1501',  'DRB', 1 )
#' predictBinding('ITKLGAKPDGKTDCT', 'DQB10201',  'DQ', 2)
#' @return  A vector cantaining four items: the peptide, the mhc-ii molecule
#' name, the probability of the binding, and the binary value (0 or 1),
#'  respectively.
#'
#' @export
predictBinding <- function(pep, mhc_name, locus,
    opt = 1) {
    locus <- toupper(locus)
    if (!(locus %in% c("DRB", "DP", "DQ"))) {
        stop("Please choose a correct locus: DRB, DP, or DQ")
    }
    v <- getInputVector(pep, mhc_name, locus)

    # w <- ifelse(opt == 1, system.file("extdata",
    #     paste(mhc_name, ".Rda", sep = ""), package = getPackageName()),
    #     system.file("extdata", paste("Q", locus,
    #         ".Rda", sep = ""), package = getPackageName()))
    w <- ifelse(opt == 1, readRDS(paste("./inst/extdata/",
        mhc_name, ".Rda", sep = "")), readRDS(paste("./inst/extdata/",
        "Q", locus, ".Rda", sep = "")))
    sumEnergy <- exp(sum(v * w))
    bindingProbability <- 1/(1 + sumEnergy)
    binaryBindingProbability <- ifelse(bindingProbability <
        0.5, 0, 1)
    return(c(pep, mhc_name, bindingProbability,
        binaryBindingProbability))
}


#' Make predictions for set of peptides and mhc-ii molecules
#'
#' This function predicts binding affinities for a dataframe  containing set
#' of peptides with set of the targeted mhcii Molcules, associated with their
#' perspective loci. It generalization of predictionBinding function.
#' @param dataset Dataframe of three columns; containing the peptide, name
#' of the mhc-ii molecule, and the name of locus, respectively. The order of
#' columns must be maintained.
#' @param opt Determines the prediction type (intra-allele or trans-alleles).
#' 1 stands for intra-alleleic and 2 for trans-allelic. The default value
#' of this is 1.
#' @return A dataframe contains four columns: the peptide, the mhc-ii molecule
#' name, the probability of the binding, and the binary value (0 or 1)
#' respectively.
#' @examples
#'  predictBindingSet <- function(dataset, opt =2)
#'
#' @export
predictBindingSet <- function(dataset, opt = 1) {
    df <- data.frame(matrix(0, nrow = dim(dataset)[1],
        ncol = 4))
    colnames(df) <- c("Peptide", "MHC Molecule",
        "Predicted Affinity", "Binary Value")
    for (i in 1:dim(dataset)[1]) {
        df[i, ] <- predictBinding(dataset[i, 1],
            dataset[i, 2], dataset[i, 3], opt)
    }
    return(df)
}
# ----------------------------------------------------------------------------
