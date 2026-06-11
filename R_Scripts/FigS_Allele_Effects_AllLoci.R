#!/usr/bin/env Rscript
# FigS_Allele_Effects_AllLoci.R
# Allele dosage effect plots for:
#   Part 1 — Chr.05 candidate gene sub-cluster peak markers
#   Part 2 — Secondary association loci (Tables S4 + S5)
#   Part 3 — yr.2025 stress-responsive loci
# Each locus: 4 panels — BLUE_exc.25 | yr.23 | yr.24 | yr.25
# yr.25 primary loci have panels reordered: yr.25 | BLUE | yr.23 | yr.24

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot); library(GWASpoly)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS")
THRESH  <- 5.52

# ── Locus definitions ─────────────────────────────────────────────────────────
# Each: list(label, marker, trait, score_blue, pve, type)
# type: "blue" = BLUE primary; "yr25" = yr.25 primary

chr05_loci <- list(
  list(label="Chr.05:47.59 Mb\n(g54676, ADK)",    marker="Chr.05_47585816", trait="DTFlower", score=10.033, pve=6.51, type="blue"),
  list(label="Chr.05:47.72 Mb\n(g54687, FIGL1)",  marker="Chr.05_47722057", trait="DTFlower", score=5.761,  pve=3.67, type="blue"),
  list(label="Chr.05:47.93 Mb\n(g54700, MFS)",    marker="Chr.05_47926790", trait="DTFlower", score=11.039, pve=7.22, type="blue"),
  list(label="Chr.05:48.11 Mb\n(g54717, ACT)",    marker="Chr.05_48114689", trait="DTFruit",  score=10.053, pve=6.52, type="blue"),
  list(label="Chr.05:48.17 Mb\n(g54726, Exo70)",  marker="Chr.05_48169153", trait="DTFlower", score=8.96,   pve=5.79, type="blue"),
  list(label="Chr.05:48.19 Mb\n(g54728, DHHC)",   marker="Chr.05_48188630", trait="DTFlower", score=9.378,  pve=6.09, type="blue"),
  list(label="Chr.05:48.24 Mb\n(g54732, BAR)",    marker="Chr.05_48243619", trait="DTFruit",  score=11.736, pve=7.69, type="blue"),
  list(label="Chr.05:48.25 Mb\n(g54732, BAR peak)*",marker="Chr.05_48252083",trait="DTFlower",score=11.86,  pve=7.79, type="blue"),
  list(label="Chr.05:48.49 Mb\n(g54753, H3K9)",   marker="Chr.05_48486332", trait="DTFlower", score=6.851,  pve=4.35, type="blue"),
  list(label="Chr.05:48.62 Mb\n(g54758, FANCM)",  marker="Chr.05_48616024", trait="DTFruit",  score=9.441,  pve=6.14, type="blue")
)

secondary_loci <- list(
  list(label="Chr.06:22–23 Mb†\n(g65194, RING-Zn)", marker="Chr.06_22292639", trait="DTFlower", score=6.988, pve=4.47, type="blue"),
  list(label="Chr.02:3–4 Mb\n(g13461, JmjC)",     marker="Chr.02_3615748",   trait="DTFlower", score=8.207, pve=5.39, type="blue"),
  list(label="Chr.02:37–38 Mb\n(g15412, A2 pep)", marker="Chr.02_37452394",  trait="DTFlower", score=7.718, pve=4.95, type="blue"),
  list(label="Chr.02:47–48 Mb\n(g16236, C2H2-Zn)",marker="Chr.02_47854803",  trait="DTFlower", score=6.100, pve=3.85, type="blue"),
  list(label="Chr.03:16–17 Mb\n(g29229, ABC tr)", marker="Chr.03_16371025",  trait="DTFruit",  score=7.721, pve=4.94, type="blue"),
  list(label="Chr.05:20–21 Mb\n(g53183, unanno)", marker="Chr.05_20189714",  trait="DTFruit",  score=7.537, pve=4.82, type="blue"),
  list(label="Chr.05:44–45 Mb\n(g54446, CDF tr)", marker="Chr.05_44735183",  trait="DTFruit",  score=7.359, pve=4.69, type="blue"),
  list(label="Chr.05:45–46 Mb\n(g54490, cNMP-bd)",marker="Chr.05_45296330",  trait="DTFruit",  score=8.370, pve=5.40, type="blue"),
  list(label="Chr.05:46–47 Mb\n(g54559, Gly hyd)",marker="Chr.05_46175030",  trait="DTFlower", score=6.874, pve=4.37, type="blue"),
  list(label="Chr.06:29–30 Mb\n(g65582, unanno)", marker="Chr.06_29469225",  trait="DTFlower", score=6.533, pve=4.19, type="blue"),
  list(label="Chr.06:42–43 Mb\n(g66181, TAR1)",   marker="Chr.06_42369200",  trait="DTFlower", score=6.778, pve=4.32, type="blue"),
  list(label="Chr.08:18–19 Mb\n(g87092, A2 pep)", marker="Chr.08_18849701",  trait="DTFlower", score=6.591, pve=4.22, type="blue"),
  list(label="Chr.08:27–28 Mb†\n(g87477, Expans)",marker="Chr.08_27127327",  trait="DTFlower", score=6.157, pve=3.89, type="blue"),
  list(label="Chr.08:28–29 Mb\n(g87581, GatB)",   marker="Chr.08_28429653",  trait="DTFlower", score=6.186, pve=4.00, type="blue"),
  list(label="Chr.08:29–30 Mb\n(g87637, MuDR TE)",marker="Chr.08_29470828",  trait="DTFlower", score=6.783, pve=4.36, type="blue"),
  list(label="Chr.09:1–2 Mb\n(g97157, Sm RNP)",   marker="Chr.09_1693527",   trait="DTFlower", score=6.205, pve=3.93, type="blue"),
  list(label="Chr.11:19–20 Mb†\n(g120148, Integ)",marker="Chr.11_19498635",  trait="DTFlower", score=5.852, pve=3.74, type="blue")
)

stress_loci <- list(
  list(label="Chr.11:45.5 Mb\n(DTFlower yr.25)", marker="Chr.11_45494107", trait="DTFlower", score=NA, pve=6.10, type="yr25"),
  list(label="Chr.07:25.0 Mb\n(DTFruit yr.25)",  marker="Chr.07_24954089", trait="DTFruit",  score=NA, pve=4.58, type="yr25")
)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading scans...\n")
load(file.path(BASE, "LINEAR/RedoNo25/GWASpoly_scans.reBLUE.RData")); scan_exc25 <- data.loco.scan
load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"));       scan_all   <- data.loco.scan



cat("Loading genotypes...\n")
load(file.path(BASE, "LINEAR/RedoNo25/data_loco.tet.reBLUE.RData"))
geno_mat <- data.loco@geno; sample_ids <- rownames(geno_mat)
rm(data.loco); gc()

cat("Loading phenotypes...\n")
pheno_fl_exc <- fread(file.path(BASE,"Phenotypes/BLUE_exc.25/DTFlower_2324only_BLUE.txt"))
pheno_fr_exc <- fread(file.path(BASE,"Phenotypes/BLUE_exc.25/DTFruit_2324only_BLUE.txt"))
setnames(pheno_fl_exc,1:2,c("Sample","Pheno"))
setnames(pheno_fr_exc,1:2,c("Sample","Pheno"))
pheno_all <- fread(file.path(BASE,"Phenotypes/SHB_BLUE_all_pheno.csv"))
setnames(pheno_all,"DNA ID","Sample")

# ── Helper functions ───────────────────────────────────────────────────────────
get_score <- function(scan, marker, trait) {
  s <- scan@scores[[trait]]
  if (is.null(s) || !marker %in% rownames(s)) return(NA_real_)
  round(max(as.numeric(s[marker,,drop=FALSE]),na.rm=TRUE),2)
}

get_dosage <- function(marker) {
  if (!marker %in% colnames(geno_mat)) return(NULL)
  dos <- as.numeric(geno_mat[,marker])
  data.table(Sample=sample_ids, Dosage=dos,
             DosClass=factor(floor(pmin(dos+0.5,4)),levels=0:4))
}

get_cld <- function(df) {
  tryCatch({
    df2 <- copy(df); df2[,DosChar:=as.character(DosClass)]
    grp_ok <- df2[,.N,by=DosChar][N>=2]
    if (nrow(grp_ok)<2) return(NULL)
    df_sub <- df2[DosChar %in% grp_ok$DosChar]; df_sub[,DosChar:=factor(DosChar)]
    mod   <- aov(Pheno~DosChar,data=df_sub)
    tukey <- TukeyHSD(mod,"DosChar")$DosChar
    groups <- levels(df_sub$DosChar); n <- length(groups)
    adj <- matrix(TRUE,n,n,dimnames=list(groups,groups))
    for (nm in names(tukey[,"p adj"])) {
      pts <- strsplit(nm,"-")[[1]]
      if (length(pts)==2&&all(pts%in%groups)){
        adj[pts[1],pts[2]] <- tukey[nm,"p adj"]>=0.05
        adj[pts[2],pts[1]] <- tukey[nm,"p adj"]>=0.05
      }
    }
    ls <- setNames(vector("list",n),groups); cur <- 1L
    for (g in groups) {
      added <- FALSE
      for (let in unique(unlist(ls))) {
        mem <- groups[sapply(groups,function(x)let %in% ls[[x]])]
        if (all(adj[g,mem])){ls[[g]]<-c(ls[[g]],let);added<-TRUE}
      }
      if (!added||length(ls[[g]])==0){ls[[g]]<-c(ls[[g]],letters[cur]);cur<-cur+1L}
    }
    cld <- sapply(groups,function(g)paste(sort(unique(ls[[g]])),collapse=""))
    wt  <- df_sub[,.(w_top=boxplot.stats(Pheno)$stats[5],
                     DosChar=as.character(DosClass[1])),by=DosClass][,.(DosChar=as.character(DosClass),w_top)]
    cld_dt <- merge(data.table(DosClass=groups,cld_letter=cld),wt,by.x="DosClass",by.y="DosChar",all.x=TRUE)
    cld_dt[,label_y:=w_top+diff(range(df$Pheno,na.rm=TRUE))*0.10]
    cld_dt
  }, error=function(e) NULL)
}

pub_theme <- theme_classic(base_size=10) +
  theme(axis.text=element_text(color="black",size=8.5),
        axis.title=element_text(size=9.5,face="bold"),
        plot.title=element_text(size=9,face="bold"),
        plot.subtitle=element_text(size=7.5,color="grey40"),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5))

COL <- list(
  DTFlower=list(blue="#8B0000", yr23="#1A6B3C", yr24="#5B2C8D", yr25="#B45309"),
  DTFruit =list(blue="#1565C0", yr23="#0D7A3C", yr24="#4A148C", yr25="#8B4513")
)

make_box <- function(dos_dt, pheno_dt, pheno_col, fill_col,
                     title_txt, sub_txt, is_primary=FALSE, y_label=NULL) {
  df <- merge(dos_dt, pheno_dt[,.(Sample,Pheno=get(pheno_col))], by="Sample")
  df <- df[!is.na(Dosage)&!is.na(Pheno)]
  if (nrow(df)<5) return(NULL)
  counts <- df[,.N,by=DosClass][order(DosClass)]
  yr <- diff(range(df$Pheno,na.rm=TRUE))
  yb <- min(df$Pheno,na.rm=TRUE)-yr*0.22
  cld_dt <- get_cld(df)
  yt <- if(!is.null(cld_dt)) max(cld_dt$label_y,na.rm=TRUE)+yr*0.05 else max(df$Pheno,na.rm=TRUE)+yr*0.1

  p <- ggplot(df,aes(x=DosClass,y=Pheno)) +
    geom_boxplot(fill=fill_col,colour="#1a1a1a",alpha=0.75,
                 outlier.size=0.4,outlier.alpha=0.3,width=0.52) +
    geom_jitter(width=0.14,alpha=0.14,size=0.4,colour="#1a1a1a") +
    stat_summary(fun=mean,geom="point",shape=23,fill="white",
                 colour=fill_col,size=2.0,stroke=1.0) +
    geom_text(data=counts,
              aes(x=DosClass,y=yb+yr*0.08,label=paste0("n=",N)),
              size=1.8,colour="grey40",angle=90,hjust=0,vjust=0.5,
              inherit.aes=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(limits=c(yb,yt)) +
    labs(title=title_txt,subtitle=sub_txt,x="ALT copies",
         y=if(!is.null(y_label))y_label else NULL) +
    pub_theme +
    theme(plot.title=element_text(size=9,face="bold",colour=fill_col),
          plot.subtitle=element_text(size=7,colour="grey40"),
          panel.border=element_rect(color=fill_col,fill=NA,
                                    linewidth=if(is_primary)1.2 else 0.6))
  if (!is.null(cld_dt)&&nrow(cld_dt)>0)
    p <- p+geom_text(data=cld_dt,aes(x=DosClass,y=label_y,label=cld_letter),
                     size=3.0,fontface="bold",colour=fill_col,inherit.aes=FALSE)
  p
}

build_row <- function(loc) {
  trait <- loc$trait
  cols  <- COL[[trait]]
  ylab  <- if(trait=="DTFlower")"Days to 50% Flowering" else "Days to 50% Ripe Fruit"
  col23 <- paste0(trait,"_yr.23"); col24 <- paste0(trait,"_yr.24"); col25 <- paste0(trait,"_yr.25")
  pheno_exc <- if(trait=="DTFlower") pheno_fl_exc else pheno_fr_exc

  dos_dt <- get_dosage(loc$marker)
  if (is.null(dos_dt)){cat("  Missing:",loc$marker,"\n"); return(NULL)}

  sc_exc <- if(!is.na(loc$score)) loc$score else get_score(scan_exc25, loc$marker, paste0(trait,"_BLUE"))
  sc_23  <- get_score(scan_all,  loc$marker, col23)
  sc_24  <- get_score(scan_all,  loc$marker, col24)
  sc_25  <- get_score(scan_all,  loc$marker, col25)

  fmt_sub <- function(sc, pve=NULL) {
    parts <- character(0)
    if(!is.na(sc)) parts <- c(parts, sprintf("-log₁₀p=%.2f",sc))
    if(!is.null(pve)&&!is.na(pve)) parts <- c(parts, sprintf("PVE=%.1f%%",pve))
    paste(parts,collapse="  |  ")
  }

  if (loc$type == "blue") {
    p1 <- make_box(dos_dt,pheno_exc,"Pheno",       cols$blue,"BLUE exc.25",fmt_sub(sc_exc,loc$pve),TRUE,ylab)
    p2 <- make_box(dos_dt,pheno_all,col23,          cols$yr23,"yr.2023",    fmt_sub(sc_23))
    p3 <- make_box(dos_dt,pheno_all,col24,          cols$yr24,"yr.2024",    fmt_sub(sc_24))
    p4 <- make_box(dos_dt,pheno_all,col25,          cols$yr25,"yr.2025",    fmt_sub(sc_25))
  } else {
    # yr.25 primary: lead with yr.25
    p1 <- make_box(dos_dt,pheno_all,col25,          cols$yr25,"yr.2025",    fmt_sub(sc_25,loc$pve),TRUE,ylab)
    p2 <- make_box(dos_dt,pheno_exc,"Pheno",        cols$blue,"BLUE exc.25",fmt_sub(sc_exc))
    p3 <- make_box(dos_dt,pheno_all,col23,          cols$yr23,"yr.2023",    fmt_sub(sc_23))
    p4 <- make_box(dos_dt,pheno_all,col24,          cols$yr24,"yr.2024",    fmt_sub(sc_24))
  }

  panels <- Filter(Negate(is.null),list(p1,p2,p3,p4))
  if (length(panels)==0) return(NULL)

  row_plot <- plot_grid(plotlist=panels, nrow=1,
                        rel_widths=c(1.12,rep(1,length(panels)-1)),
                        align="h",axis="tb")

  lbl <- ggdraw()+draw_label(loc$label,fontface="bold",size=7.5,angle=90,x=0.5,y=0.5)
  plot_grid(lbl, row_plot, nrow=1, rel_widths=c(0.07,1))
}

# ── Build all rows ─────────────────────────────────────────────────────────────
cat("Building panels...\n")
save_section <- function(locus_list, section_name, n_per_page=4) {
  rows <- lapply(locus_list, function(loc){
    cat(sprintf("  [%s] %s\n", section_name, gsub("\n"," ",loc$label)))
    build_row(loc)
  })
  rows <- Filter(Negate(is.null), rows)
  split(rows, ceiling(seq_along(rows)/n_per_page))
}

pages_chr05 <- save_section(chr05_loci,     "Chr05")
pages_sec   <- save_section(secondary_loci, "Secondary")
pages_str   <- save_section(stress_loci,    "Stress", n_per_page=2)

# ── Save PDFs ─────────────────────────────────────────────────────────────────
save_pages <- function(pages, pdf_path, png_path, w=11, h=12) {
  pdf(pdf_path, width=w, height=h, onefile=TRUE)
  for (pg in pages) {
    rh <- rep(1,length(pg))
    print(plot_grid(plotlist=pg,ncol=1,rel_heights=rh))
  }
  dev.off()
  ggsave(png_path, plot_grid(plotlist=pages[[1]],ncol=1,rel_heights=rep(1,length(pages[[1]]))),
         width=w, height=h, dpi=300)
  cat("Saved:",pdf_path,"\n")
}

save_pages(pages_chr05,
  file.path(OUT_DIR,"FigS_Allele_Effects_Chr05_Subclusters.pdf"),
  file.path(OUT_DIR,"FigS_Allele_Effects_Chr05_Subclusters_p1.png"))

save_pages(pages_sec,
  file.path(OUT_DIR,"FigS_Allele_Effects_Secondary_Loci.pdf"),
  file.path(OUT_DIR,"FigS_Allele_Effects_Secondary_Loci_p1.png"))

save_pages(pages_str,
  file.path(OUT_DIR,"FigS_Allele_Effects_Stress_Loci.pdf"),
  file.path(OUT_DIR,"FigS_Allele_Effects_Stress_Loci_p1.png"),
  w=11, h=7)

cat("\nAll done.\n")
