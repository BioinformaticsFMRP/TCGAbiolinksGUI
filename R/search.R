
# internal: used by biOmicsSearch
systemSearch <- function(term,env) {
    success <- get("success", envir = env)
    ont <- get("ont", envir = env)
    systems <- get("systems", envir = as.environment("package:biOmics"))

    if (!success) {
        found <- intersect(term, systems$BTO)
        if (length(found) > 0) {
            assign("success", TRUE, envir = env)
            assign("solution", found, envir = env)
            return()
        } else {
            parentTerm <- rols::parents(term[1], ont)
            parentTerm <- names(parentTerm)
            sapply(parentTerm, function(x) {
                systemSearch(x,env)
            })
            return()
        }
    }
}

map.to.bto <- function(term,env) {
    success <- get("success", envir = env)
    # Does the term exists in BTO?
    aux <- term
    names(aux) <- NULL
    query <- rols::olsQuery(aux, "BTO", exact = FALSE)

    if (length(query) > 0) {
        # term found in BTO
        term.found <- names(query)
        sapply(term.found, function(x) {
            if (!success) {
                systemSearch(x, env)
            }
        })
    }
    if (success) {
        return()
    } else {
        ontology <- unlist(strsplit(names(term), ":"))[1]
        parentTerm <- rols::parents(names(term), ontology)

        for (i in seq(parentTerm)) {
            map.to.bto(parentTerm[i])
        }
    }
}


is.mapped <- function(term,env) {

    # search in roamap Could be - Accession, Sample.Name,
    # Experiment
    idx <- grep(term, biosample.roadmap$biosample, ignore.case = TRUE)
    if (length(idx) > 0) {
        message("found in roadmap")
        res <- biosample.roadmap[idx[1], ]$BTO
        if (!is.na(res)) {
            assign('success', TRUE, envir = env)
            assign('solution', res, envir = env)
            return()
        }
    }
    # search in encode
    idx <- grep(term, biosample.encode$biosample, ignore.case = TRUE)

    if (length(idx) > 0) {
        message("found in encode")
        res <- biosample.encode[idx[1], ]$BTO
        if (!is.na(res)) {
            assign('success', TRUE, envir = env)
            assign('solution', res, envir = env)

            return()
        }

    }
    # search in tcga
    idx <- grep(term, biosample.tcga$biosample, fixed = TRUE)
    if (length(idx) == 0) {
        idx <- grep(term, biosample.tcga$biosample, ignore.case = TRUE)
    }
    if (length(idx) > 0) {
        message("found in TCGA")
        res <- biosample.tcga[idx[1], ]$BTO
        assign('success', TRUE, envir = env)
        assign('solution', res, envir = env)
        return()
    }
}

#' @title Searches a term in TCGA, ENCODE, ROADMAP
#' @description
#'  Search a term using an ontology in the TCGA, ENCODE, ROADMAP databases
#' @param term Term to be searched. Example: 'brain' 'u87' etc.
#' @param experiment Experiment type
#' @param report Create a summary plot of the result found? Deafult: TRUE
#' @param dir.report Directory to save the html report. Default: "searchSummary"
#' \tabular{llll}{
#'Microarray \tab MiRNAMicroArray \tab RRBS \tab DNAsequencing\cr
#'ExpressionArray \tab Firehose \tab ChipSeq \tab fiveC \cr
#'ExonArray \tab DNAMethylation  \tab MRESeq  \tab RepliSeq \cr
#'RNASeq \tab miRNASeq \tab Rampage \tab Others
#'}
#' @examples inst/examples/biomicsSearch.R
#' @importFrom rols olsQuery term parents isIdObsolete
#' @importFrom TCGAbiolinks TCGAquery
#' @export
#' @return A dataframe with the results of the query if it
#'         was successful
biOmicsSearch <- function(term,
                          experiment = NULL,
                          report = FALSE,
                          dir.report = "searchSummary") {

    message(paste("biOmics is searching for:", term, "\nSearching..."))
    start.time <- Sys.time()
    env <- new.env()
    assign('success', FALSE, envir = env)
    success <- get("success", envir = env)
    assign('solution',FALSE, envir = env )
    assign('exper', experiment, envir = env)

    # Step 0: verify if term is valid.
    if (!is.valid.term(term)) {
        return()
    }

    # Step 0.5: verify if experiment is valid.
    if (!is.experiment(experiment)) {
        return()
    }

    # Step 1: verify if term search has been mapped by us.
    if (!success) {
        is.mapped(term, env)
        success <- get("success", envir = env)
        solution <- get("solution", envir = env)
    }

    # Step 2: search for the term in the BTO ontology.
    if (!success) {
        message("Not found in the cache, searching in the ontology...")
        # Term has not been mapped before.
        assign("ont", "BTO" , envir = env)
        ont <- get("ont", envir = env)
        query <- rols::olsQuery(term, ont, exact = TRUE)

        if (length(query) > 0) {
            term.found <- names(query)
            systemSearch(term.found,env)
        } else {
            query <- rols::olsQuery(term, ont, exact = FALSE)
            if (length(query) > 0) {
                term.found <- names(query)
                sapply(term.found, function(x) {
                    if (!rols::isIdObsolete(x, ont)) {
                        systemSearch(x,env)
                    }
                })
            }
        }

        # Not found in bto or cache we are going to search in other
        # ontologies and see if from the other ontologies we can
        # reach BTO
        if (!success) {
            message("Not found in the cache, not found in BTO...")

            ont.vec <- c("EFO", "CL", "UBERON")
            for (i in ont.vec) {
                query <- rols::olsQuery(term, i)
                if (length(query) > 0) {
                    # term was found in other ontology!
                    for (j in seq(query)) {
                        map.to.bto(query[j],env)
                        success <- get("success", envir = env)
                        if (success)  break
                    }
                    if (success) {
                        message(paste0("Found in:",i,"!"))
                        break
                    }
                }
            }
        }
        success <- get("success", envir = env)
        solution <- get("solution", envir = env)
    }

    end.time <- Sys.time()
    time.taken <- end.time - start.time
    message(paste("Time taken: ", round(time.taken, 2), "s \n"))

    if (success) {
        return(showResults(solution, experiment, report, dir.report))
    }

    message("Term not found")
    return(NULL)
}

# show the results to the user
#' @import ggplot2
#' @keywords internal
showResults <- function(solution, exper, report = FALSE, path) {

    # Get the samples that matches the result of the query
    # Databases were matched manually to systems
    pat <- unlist(strsplit(solution[1], ","))

    idx <- apply(sapply(pat, function(x) {
        grepl(x, biosample.encode$BTO)
    }), 1, any)
    enc.samples <- biosample.encode[idx, ]$biosample
    idx <- apply(sapply(pat, function(x) {
        grepl(x, biosample.roadmap$BTO)
    }), 1, any)
    rmap.samples <- biosample.roadmap[idx, ]$biosample
    idx <- apply(sapply(pat, function(x) {
        grepl(x, biosample.tcga$BTO)
    }), 1, any)
    tcga.samples <- biosample.tcga[idx, ]$biosample

    # Select the database rows
    enc.result <- encode.db[is.element(encode.db$biosample, enc.samples),]
    rmap.result <- roadmap.db[is.element(roadmap.db$`Standardised epigenome name`,
                                         rmap.samples), ]
    disease <- sapply(strsplit(tcga.samples, split = " - "),
                      function(x) {x[1]})
    tcga.result <- TCGAquery(tumor = disease, level = 3)

    # Select experiments
    if (! is.null(exper)) {
        message("Filtering by experiment")
        idx <- apply(sapply(platforms[grep(exper, platforms$Standard,
                                           ignore.case = TRUE), 2],
                            function(x) {
                                grepl(x,
                                      enc.result$assay,
                                      ignore.case = TRUE)
                            }
        ), 1, any)
        enc.result <- enc.result[idx, ]

        idx <- apply(sapply(platforms[grep(exper, platforms$Standard,
                                           ignore.case = TRUE), 2],
                            function(x) {
                                grepl(x,
                                      rmap.result$MARK,
                                      ignore.case = TRUE)
                            }
        ), 1, any)
        rmap.result <- rmap.result[idx, ]

        idx <- apply(sapply(platforms[grep(exper, platforms$Standard,
                                           ignore.case = TRUE), 2],
                            function(x) {
                                grepl(x,
                                      tcga.result$baseName,
                                      ignore.case = TRUE)
                            }
        ), 1, any)
        if (length(idx) > 0) {
            tcga.result <- tcga.result[idx, ]
        } else {
            tcga.result <- tcga.result[1:nrow(tcga.result), ]
        }

    }

    # Preparing the output table
    colnames(rmap.result)[c(2,4,1)] <- c("ID", "Sample", "Experiment")
    colnames(enc.result)[1:3] <- c("ID", "Sample", "Experiment")
    colnames(tcga.result)[c(7, 11, 10)] <- c("ID", "Sample", "Experiment")
    rmap.result <- as.data.frame(rmap.result)
    results <- NULL
    if(nrow(enc.result) > 0 & is.null(results))  {
        results <- enc.result[,1:3]
        if(nrow(rmap.result) > 0 & !is.null(results)) results <- rbind(results, rmap.result[,c(2,4,1)])
        if(nrow(tcga.result) > 0 & !is.null(results)) results <- rbind(results, tcga.result[,c(7,11, 10)])
    } else if(nrow(rmap.result) > 0 & is.null(results)) {
        results <- rmap.result
        if(nrow(tcga.result) > 0 & !is.null(results)) results <- rbind(results, tcga.result[,c(7,11, 10)])
    } else if(nrow(tcga.result) > 0 & is.null(results)) {
        results <- tcga.result
    }


    database <- c(rep("encode", nrow(enc.result)),
                  rep("roadmap", nrow(rmap.result)),
                  rep("tcga", nrow(tcga.result)))
    results <- cbind(database, results)

    if (report) {
        create.report(results,path = path,
                      system = sapply(pat, function(x) {subset(systems, systems$BTO == x)$system}))
    }
    return(results)
}

create.summary.plot <- function(data, col, title, filename, path){
    .e <- environment()

    g <- ggplot(data, aes(factor(database),
                          fill = data[,col]),
                environment = .e) +
        geom_bar(position = "fill") +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.key = element_rect(colour = 'white'),
              legend.justification=c(1,1), legend.position=c(1,1)) +
        ggtitle(title) +
        labs(x="Database", y="Percentage of samples") +
        scale_fill_discrete(name=col) +
        theme(legend.direction ="vertical",legend.position = "bottom")+
        guides(fill=guide_legend(ncol=4))

    plot(g)
    ggsave(g, filename = file.path(path, filename),
           height = 5, width = 5, scale = 1.5)

}

#' @importFrom knitr kable
is.experiment <- function(experiment) {

    v <- unique(platforms$Standard)
    if(is.null(experiment)){
        return(TRUE)
    }
    if ((length(grep(experiment, v, ignore.case = TRUE)) > 0) &
        (nchar(experiment) >= 3)) {
        return(TRUE)
    } else {
        df <- as.data.frame(matrix(sort(unique(platforms$Standard)),
                                   ncol = 3))
        print(kable(df, col.names = NULL, format = "pandoc",
                    caption = "Experiment"))
        cat("=======================================================\n")
        cat("ERROR: Experiment not found. Select from the table above.\n")
        cat("       Select from the table above.\n")
        cat("=======================================================\n")
        return(FALSE)
    }
}

is.valid.term <- function(term) {
    if (nchar(term) >= 3) {
        return(TRUE)
    } else {
        message(paste0("ERROR: ", term, " is not valid.",
                       " Specify a term of at least 3 characters"))
        return(FALSE)
    }
}

# ------------------------ Roadmap search
#' @title Roadmap search
#' @description
#'    Searches in the roadmap database
#'    Source: https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=14
#' @param sample Example:
#' \tabular{ll}{
#'H1 cell line  \tab Fetal Heart\cr
#'IMR90 fetal lung fibroblasts Cell Line \tab Fetal Kidney\cr
#'Aorta \tab lung\cr
#'Fetal Adrenal Gland \tab CD34 cultured cells
#' }
#' @param experiment Examples:
#'\tabular{llll}{
#'DNase   \tab H2A.Z  \tab H2AK5ac \tab H2BK120ac   \cr
#'H2BK12ac   \tab H2BK15ac \tab H2BK20ac  \tab H2BK5ac \cr
#'H3K14ac  \tab H3K18ac   \tab H3K23ac  \tab H3K23me2  \cr
#'H3K27ac    \tab H3K27me3  \tab H3K36me3   \tab  H3K4ac   \cr
#'H3K4me1   \tab H3K4me2   \tab H3K4me3  \tab H4K5ac
#'}
#' @export
#' @return Dataframe with the results of the query
roadmapSearch <- function(sample = NULL,
                          experiment = NULL) {

    db <- roadmap.db

    if (!is.null(sample)) {
        id <- sapply(sample, function(x) {
            grepl(x,db$Standardised.epigenome.name,ignore.case = T)
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(experiment)) {
        id <- sapply(experiment, function(x) {
            db$MARK == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    return(db)
}


# ------------------------ Encode search
#' @title Encode search
#' @description
#'    Searches in the encode database
#' @param accession GEO sample accession Ex: ENCSR337VPD
#' @param biosample Example:
#'\tabular{lllll}{
#'  hepatocyte             \tab GM12878    \tab induced pluripotent stem cell
#'   \tab kidney    \tab Daoy                    \cr
#'  neural progenitor cell \tab A673       \tab bipolar spindle neuron
#'   \tab heart     \tab A172                    \cr
#'  K562                   \tab Karpas-422 \tab adipose tissue
#'  \tab liver     \tab LHCN-M2                 \cr
#'  A549                   \tab MM.1S      \tab adrenal gland
#'  \tab testis    \tab skeletal muscle myoblast\cr
#'  brain                  \tab HT1080     \tab female gonad
#'  \tab RPMI-7951 \tab myotube                 \cr
#'  spleen                 \tab SK-N-DZ    \tab lung
#'  \tab M059J     \tab H7-hESC                 \cr
#'  pancreas               \tab SK-MEL-5   \tab sigmoid colon
#'  \tab H4        \tab astrocyte               \cr
#'  HepG2                  \tab NCI-H460   \tab small intestine
#'  \tab SJSA1     \tab HeLa-S3
#'}
#' @param assay Example:
#'\tabular{ll}{
#'  RAMPAGE  \tab RRBS \cr
#'  ChIP-seq \tab FAIRE-seq \cr
#'  RNA-seq  \tab Repli-seq \cr
#'  shRNA knockdown followed by RNA-seq
#'  \tab DNA methylation profiling by array assay            \cr
#'  DNase-seq \tab RIP-seq \cr
#'  iCLIP     \tab Switchgear \cr
#'  ChIA-PET  \tab CAGE   \cr
#'  MeDIP-seq assay \tab RIP-chip   \cr
#'  MRE-seq         \tab Repli-chip \cr
#'  RNA profiling by array assay
#'   \tab whole-genome shotgun bisulfite sequencing           \cr
#'  RNA-PET         \tab MNase-seq  \cr
#'  5C              \tab DNA-PET    \cr
#'  comparative genomic hybridization by array
#'  \tab protein sequencing by tandem mass spectrometry assay
#'}
#' @param lab Example:
#' \tabular{lllll}{
#'Thomas Gingeras, CSHL    \tab John Stamatoyannopoulos, UW
#' \tab Kevin Struhl, HMS    \tab Gregory Crawford, Duke
#'  \tab David Gilbert, FSU      \cr
#'Kevin White, UChicago    \tab Barbara Wold, Caltech
#'\tab Peggy Farnham, USC   \tab Morgan Giddings, UNC
#'\tab Ali Mortazavi, UCI      \cr
#'Michael Snyder, Stanford \tab Gene Yeo, UCSD
#'\tab Job Dekker, UMass    \tab Scott Tenenbaum, SUNY-Albany
#'\tab Ross Hardison, PennState\cr
#'Brenton Graveley, UConn  \tab Bradley Bernstein, Broad
#'\tab Vishwanath Iyer, UTA \tab Piero Carninci, RIKEN
#'\tab Bing Ren, UCSD          \cr
#'Richard Myers, HAIB      \tab Yijun Ruan, GIS
#'\tab Jason Lieb, UNC      \tab Sherman Weissman, Yale \tab Joe Ecker, Salk
#'}
#' @param target Target Example:
#'\tabular{lllll}{
#'  Control            \tab SERBP1 \tab FAM120A \tab DDX27  \tab RPLP0 \cr
#'  rabbit-IgG-control \tab RRP9   \tab EIF3D   \tab CSTF2  \tab RPL23A\cr
#'  NUSAP1             \tab RPS19  \tab EFTUD2  \tab BCLAF1 \tab RCC2  \cr
#'  XRN2               \tab PRPF8  \tab EEF2    \tab BCCIP  \tab RBM27 \cr
#'  UCHL5              \tab PES1   \tab DDX55   \tab UPF2   \tab PSIP1 \cr
#'  TFIP11             \tab PA2G4  \tab DDX52   \tab TROVE2 \tab PKM   \cr
#'  SUPV3L1            \tab NONO   \tab DDX51   \tab RPS5   \tab PHF6  \cr
#'  SRP68              \tab GEMIN5 \tab DDX28   \tab RPS2   \tab
#'  Non-specific target control
#'}
#' @param description Description of the sample
#' @param organism The organism that the sample belongs to
#'\tabular{l}{
#'Homo sapiens \cr Mus musculus \cr Drosophila melanogaster
#'}
#' @example inst/examples/encodeSearch.R
#' @import ggplot2
#' @export
#' @return Dataframe with the query result
encodeSearch <- function(accession = NULL,
                         biosample = NULL,
                         assay = NULL,
                         lab = NULL,
                         target = NULL,
                         description = NULL,
                         organism = NULL,
                         exact = TRUE) {

    # roadmap.verify.input(GEO.Accession,sample,experiment,
    # center,embargo.end.date)
    valid <- validadeEncode(accession, biosample, assay, lab, target,
                            description, organism)
    if(!(valid)){
        return(NULL)
    }

    db <- encode.db

    if (!is.null(biosample)) {
        id <- sapply(biosample, function(x) {
            db$biosample == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(assay)) {
        id <- sapply(assay, function(x) {
            db$assay == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(accession)) {
        id <- sapply(accession, function(x) {
            db$accession == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(lab)) {
        id <- sapply(lab, function(x) {
            db$lab == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(target)) {
        id <- sapply(target, function(x) {
            db$target == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(organism)) {
        id <- sapply(organism, function(x) {
            db$organism == x
        })
        id <- apply(id, 1, any)
        db <- db[id, ]
    }
    if (!is.null(description)) {
        id <- grep(description, db$description, ignore.case = FALSE)
        db <- db[id, ]
    }

    return(db)
}

validadeEncode <- function(accession = NULL, biosample = NULL,
                           assay = NULL, lab = NULL, target = NULL,
                           description = NULL, organism = NULL
){

    db <- encode.db

    if (!is.null(accession)) {
        if (!length(grep(accession, db$accession,
                         ignore.case = TRUE)) > 0 ){
            cat("=======================================================\n")
            cat("ERROR: Acession not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(FALSE)
        }
    }

    if (!is.null(biosample)) {
        if (!length(grep(accession, db$biosample,
                         ignore.case = TRUE)) > 0 ){
            df <- as.data.frame(matrix(sort(unique(db$biosample)),
                                       ncol = 3))
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "Encode samples"))
            cat("=======================================================\n")
            cat("ERROR: Samples not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(FALSE)
        }
    }

    if (!is.null(assay)) {
        if  (!length(grep(assay, db$assay,
                          ignore.case = TRUE)) > 0 ){
            df <- as.data.frame(matrix(sort(unique(db$assay)),
                                       ncol = 3))
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "Encode assay"))
            cat("==========================================================\n")
            cat("ERROR: Assay not found. Select from the table above.\n")
            cat("==========================================================\n")
            return(FALSE)
        }
    }

    if (!is.null(lab)) {
        if  (!length(grep(center, db$lab, ignore.case = TRUE)) > 0 ){
            df <- as.data.frame(matrix(sort(unique(db$lab)),
                                       ncol = 3))
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "Encode lab"))
            cat("==========================================================\n")
            cat("ERROR: Lab not found. Select from the table above.\n")
            cat("==========================================================\n")
            return(FALSE)
        }
    }

    if (!is.null(target)){
        if (!length(grep(target, db$target,
                         ignore.case = TRUE)) > 0 ){
            df <- as.data.frame(matrix(sort(unique(db$target)),
                                       ncol = 3))
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "Encode lab"))
            cat("=======================================================\n")
            cat("ERROR: Target not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(FALSE)
        }
    }

    if (!is.null(organism)) {
        if (!length(grep(organism, db$organism,
                         ignore.case = TRUE)) > 0 ){
            df <- as.data.frame(matrix(sort(unique(db$organism)),
                                       ncol = 3))
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "Encode lab"))
            cat("=======================================================\n")
            cat("ERROR: Organism not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(FALSE)
        }
    }
    return(TRUE)
}
