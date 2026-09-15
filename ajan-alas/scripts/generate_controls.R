suppressPackageStartupMessages({library(readr); library(dplyr)})
cmd_args <- commandArgs(trailingOnly=FALSE); file_arg <- grep("^--file=",cmd_args,value=TRUE)
script_path <- normalizePath(sub("^--file=","",file_arg)); setwd(dirname(dirname(script_path)))
param_file <- "parameters/ajan-alas-m3-prior2-s20-p4-means.csv"
model_file <- "model/ajan_alas_species_model.R"
if (!file.exists(param_file)) stop("Missing ",param_file)
if (!file.exists(model_file)) stop("Missing ",model_file,". Species-level representation must be defined before simulation.")
source(model_file)
required <- c("population_order","comparison_specs","make_species_tree","migration_pairs")
miss <- required[!vapply(required,exists,logical(1),inherits=FALSE)]; if(length(miss)) stop("Model file missing: ",paste(miss,collapse=", "))
p <- read_csv(param_file,show_col_types=FALSE); dir.create("controls",showWarnings=FALSE); dir.create("imap",showWarnings=FALSE); dir.create("output/trees",recursive=TRUE,showWarnings=FALSE)
write_imap <- function(spec){f<-file.path("imap",paste0(spec$comparison,"_",spec$config,".imap.txt")); x<-if(spec$config=="aab") c(paste0("a1 ",spec$pop_a),paste0("a2 ",spec$pop_a),paste0("b1 ",spec$pop_b)) else c(paste0("a1 ",spec$pop_a),paste0("b1 ",spec$pop_b),paste0("b2 ",spec$pop_b)); writeLines(x,f); f}
make_ctl <- function(row,spec){counts<-setNames(rep(0L,length(population_order)),population_order); if(spec$config=="aab"){counts[spec$pop_a]<-2;counts[spec$pop_b]<-1}else{counts[spec$pop_a]<-1;counts[spec$pop_b]<-2}; imap<-write_imap(spec); tree<-make_species_tree(row); mig<-vapply(migration_pairs,function(pair){nm<-paste0("W_",pair[1],"_to_",pair[2]); if(!nm%in%names(row))stop("Missing parameter: ",nm); sprintf("            %s %s %.10g",pair[1],pair[2],row[[nm]])},character(1)); paste(c("seed = -1",paste("treefile =",file.path("output","trees",paste0(row$chromosome,"_",spec$comparison,"_",spec$config,".tree.txt"))),paste("Imapfile =",imap),paste("species&tree =",length(population_order),paste(population_order,collapse=" ")),paste("                ",paste(unname(counts),collapse=" ")),tree,"loci&length = 1000000 50",paste("migration =",length(migration_pairs)),mig),collapse="\n")}
for(i in seq_len(nrow(p))){row<-p[i,];for(j in seq_len(nrow(comparison_specs))){spec<-comparison_specs[j,];writeLines(make_ctl(row,spec),file.path("controls",paste0(row$chromosome,"_",spec$comparison,"_",spec$config,".ctl")))}}
message("Generated ajan-alas species-level controls")
