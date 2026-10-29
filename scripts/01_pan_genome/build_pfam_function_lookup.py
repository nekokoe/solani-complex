#!/usr/bin/env python3
"""Build a Pfam -> functional-class lookup table for this project.

The class scheme below is a summary of the functional groupings Pfam domains cover, made
once and stored on disk so that every downstream count uses the same vocabulary. The
assignment layer has to cope with three kinds of input:

  * a family in the Pfam 38.0 clan file — 25,545 of them, and every one carries a clan;
    the class is inherited from the clan, with family-level overrides where a clan pools
    domains that do different jobs (NADP_Rossmann carries both dehydrogenases and
    epimerases, TIM_barrel both glycolytic enzymes and cellulases, and so on);
  * a family whose clan is uninformative but whose own description names the function;
  * a family that is genuinely uncharacterised ("Protein of unknown function (DUF…)"),
    which is kept as its own class rather than forced into a functional one.

Every row records how it was classified, so the table can be audited (`method` column:
clan / clan_override / keyword / unknown / unmatched).

Inputs
    /share/pasteur/database/Pfam/Pfam38.0/pfamA.txt.gz          accession, name, description
    /share/pasteur/database/Pfam/Pfam38.0/Pfam-A.clans.tsv.gz   accession, clan id/name, name, description
Outputs
    02.annotation/pfam_function_lookup.tsv
    02.annotation/pfam_function_classes.tsv
"""

import csv
import gzip
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PFAM_DIR = "/share/pasteur/database/Pfam/Pfam38.0"
PfamA = os.path.join(PFAM_DIR, "pfamA.txt.gz")
CLANS = os.path.join(PFAM_DIR, "Pfam-A.clans.tsv.gz")
OUT = os.path.join(ROOT, "02.annotation/pfam_function_lookup.tsv")
OUT_CLASSES = os.path.join(ROOT, "02.annotation/pfam_function_classes.tsv")

# ---------------------------------------------------------------- class scheme
# Each entry: class -> what belongs there. Kept short on purpose: the point is a stable
# vocabulary for counting, not an ontology.
CLASSES = [
    ("Transport", "transporters, permeases, channels, solute-binding proteins"),
    ("Transcription and chromatin", "RNA polymerase, transcription factors, histones, chromatin modifiers"),
    ("Translation and ribosome biogenesis", "ribosomal proteins, tRNA synthetases, translation factors"),
    ("RNA processing and RNA binding", "RNA recognition motifs, RNA helicases, splicing, modification"),
    ("DNA replication, repair and recombination", "polymerases, helicases, topoisomerases, transposases, ligases"),
    ("Protein folding, modification and turnover", "chaperones, ubiquitin system, proteasome, peptidyl-prolyl isomerases"),
    ("Cytoskeleton and motility", "actin, tubulin, motors, profilin, septins"),
    ("Signal transduction", "protein kinases, phosphatases with regulatory roles, G proteins, second-messenger domains"),
    ("Oxidoreduction", "dehydrogenases, oxidases, oxygenases, reductases, peroxidases, electron carriers"),
    ("Transfer of groups", "kinases, aminotransferases, acetyl/methyl transferases, phosphorylases, synthases, ligases"),
    ("Hydrolysis", "proteases, esterases, lipases, phosphatases, nucleases, glycosidases outside CAZy curation"),
    ("Carbohydrate metabolism and CAZymes", "glycoside hydrolases, glycosyltransferases, carbohydrate esterases/lyases, sugar epimerases"),
    ("Lipid and secondary metabolism", "fatty-acid and sterol enzymes, PKS/NRPS, terpene and alkaloid enzymes"),
    ("Cell wall and membrane biogenesis", "chitin/glucan synthases, GPI anchor, mannosyltransferases"),
    ("Trafficking and secretion", "SNAREs, coats, Rab/NTF2, vacuolar sorting"),
    ("Structural repeat / protein interaction", "ankyrin, TPR, WD40, LRR, Ig-like, lectin folds"),
    ("Other enzymatic", "enzymatic folds whose substrates span classes (cupin, glyoxalase, vWA-like)"),
    ("Unknown function", "DUF families and uncharacterised proteins"),
    ("No assignment", "accession absent from the Pfam 38.0 clan file, or description too generic"),
]


# ---------------------------------------------------------------- GO MF / EC axis
# Pfam families are protein domains, so what they carry natively is a MOLECULAR FUNCTION
# (GO MF), not a pathway. This axis therefore uses the EC top classes plus the non-enzyme
# GO MF children, and is what the manuscript should quote for "what these genes do".
# KEGG level-2 stays as its own vocabulary (a pathway axis) and is not mixed into this.
MF_CLASSES = [
    ("Oxidoreductase activity (EC 1)", "redox enzymes and electron carriers"),
    ("Transferase activity (EC 2)", "kinases, aminotransferases, acetyl/methyl/glycosyl transferases, ubiquitin ligases"),
    ("Hydrolase activity (EC 3)", "proteases, esterases, lipases, phosphatases, nucleases, glycosidases, ATPases"),
    ("Lyase activity (EC 4)", "decarboxylases, aldolases, dehydratases, pectate lyases"),
    ("Isomerase activity (EC 5)", "epimerases, isomerases, mutarotases"),
    ("Ligase activity (EC 6)", "tRNA synthetases, DNA ligases, AMP-binding/ATP-grasp enzymes"),
    ("Translocase / transporter activity (EC 7)", "transporters, permeases, channels, pumps, solute-binding system components"),
    ("Transcription regulator activity", "DNA-binding transcription factors and their cofactors"),
    ("Translation regulator activity", "translation factors and ribosome-associated regulators"),
    ("Molecular function regulator", "GTPase regulators, kinase regulators, chaperone cofactors"),
    ("Chaperone activity", "protein-folding chaperones"),
    ("Structural molecule activity", "cytoskeletal and structural proteins"),
    ("Motor activity", "actin- and tubulin-based motors"),
    ("Antioxidant activity", "thioredoxin, glutaredoxin, peroxidase-like redox protection"),
    ("Binding: nucleic acid", "DNA/RNA-binding domains without catalytic assignment"),
    ("Binding: protein", "protein-interaction domains and scaffolds"),
    ("Binding: small molecule or ligand", "cofactor, lipid, calcium, solute-binding domains"),
    ("Binding: carbohydrate", "carbohydrate-binding modules and lectin domains"),
    ("Signal sensor / transducer activity", "receptors and signalling modules"),
    ("Catalytic activity, unspecified", "catalytic folds whose substrate is not assignable from the domain alone"),
    ("Unknown function", "DUF families and uncharacterised proteins"),
    ("No assignment", "family absent from Pfam 38.0 metadata"),
]

# From the semantic (KEGG-style) class already assigned, to the GO MF axis. Families that
# need an EC-level nuance are overridden below by name pattern.
SEMANTIC_TO_MF = {
    "Oxidoreduction": "Oxidoreductase activity (EC 1)",
    "Transfer of groups": "Transferase activity (EC 2)",
    "Hydrolysis": "Hydrolase activity (EC 3)",
    "Carbohydrate metabolism and CAZymes": "Transferase activity (EC 2)",
    "Lipid and secondary metabolism": "Catalytic activity, unspecified",
    "Structural repeat / protein interaction": "Binding: protein",
    "Transcription and chromatin": "Transcription regulator activity",
    "Translation and ribosome biogenesis": "Binding: nucleic acid",
    "RNA processing and RNA binding": "Binding: nucleic acid",
    "DNA replication, repair and recombination": "Hydrolase activity (EC 3)",
    "Protein folding, modification and turnover": "Chaperone activity",
    "Cytoskeleton and motility": "Structural molecule activity",
    "Signal transduction": "Signal sensor / transducer activity",
    "Transport": "Translocase / transporter activity (EC 7)",
    "Trafficking and secretion": "Binding: protein",
    "Cell wall and membrane biogenesis": "Transferase activity (EC 2)",
    "Other enzymatic": "Catalytic activity, unspecified",
    "Unknown function": "Unknown function",
    "No assignment": "No assignment",
}

# EC-level refinements that the semantic class alone cannot express.
import re as _re

# EC second level. The one-digit EC class alone merges things that behave very differently
# across the pan-genome: EC 3 holds both proteases and glycoside hydrolases, and it also
# swallows the AAA+/GTPase families, which are ATPases (EC 3.6) rather than hydrolases of a
# bond. Splitting to EC x.y is the standard EC hierarchy and keeps the axis comparable to a
# GO molecular-function annotation.
EC_PATTERNS = [
    # EC 3.6 — acting on acid anhydrides (ATPases / GTPases)
    (r"^(AAA|AAA_\d+|AAA_ ATPase|MMR_HSR1|SRP54|GTP_EFTU|GTP_EFTU_D\d|Ras|Ras_\d|"
     r"Kinesin|Myosin_tail|Dynamin_N|MCM|MCM_OB|DEAD|Helicase_C|ResIII|SNF2_N|"
     r"ABC_tran|ABC2_membrane|ABC_membrane|P-loop_NTPase|NACHT|HSP70|HATPase_c)", "EC 3.6 - ATPase / GTPase"),
    # EC 3.4 — peptidases (protein turnover)
    (r"(Peptidase|peptidase|Protease|protease|Trypsin|Trypsin_2|Lysyl-oxidase|"
     r"OTCace|Aminopeptidase|Carboxypeptidase|Serpin|Inhibitor_I\d+|"
     r"Lon_C|Clp_N|ClpB_D2-small|PUA|Propeptide_multi-cystatin)", "EC 3.4 - peptidase (protein turnover)"),
    # EC 3.2 — glycosylases (CAZymes)
    (r"(Glyco_hydro|Glycosyl_hydrolase|Glucosid|Galactosid|Mannosid|Chitin|Xylanase|"
     r"Cellulase|Arabinofuranosid|Rhamnosid|Fucosid|Sialidase|Lysozyme|Pectate|"
     r"GH\d+|CBM\d+|Glyco_hydro_\d+)", "EC 3.2 - glycosylase (carbohydrate)"),
    # EC 3.1 — esterases / lipases / phosphatases / nucleases
    (r"(Abhydrolase|COesterase|Esterase|esterase|Lipase|lipase|Thioesterase|"
     r"Phosphatase|phosphatase|Phosphoesterase|Nuclease|RNase|DNase|"
     r"Metallophosphoesterase|PAP2|Calcineurin_like|DSPc|PTPc)", "EC 3.1 - esterase / phosphatase / nuclease"),
    # EC 1.14 / 1.13 — oxygenases (P450, secondary metabolism)
    (r"^(p450|P450|Cyt_P450|FA_hydroxylase|PhyH|TauD|TauE|DAO|FMO-like|"
     r"Lycopene_cycl|Beta_elim_lyase|IPNS|Dioxygenase_N|Ring_hydroxyl_A)", "EC 1.13/1.14 - oxygenase (P450-type)"),
    # EC 1.1 — dehydrogenases acting on CH-OH
    (r"^(adh_short|adh_short_C2|ADH_zinc_N|ADH_N|Aldedh|KR|3Beta_HSD|"
     r"NAD_binding_\d+|NADP_Rossmann|GDP_Man_Dehyd|GFO_IDH_MocA|Aldo_ket_red|"
     r"SDR|Epimerase|TrkA_N)", "EC 1.1 - dehydrogenase (CH-OH)"),
    # EC 2.7 — phosphotransferases
    (r"^(Pkinase|Pkinase_Tyr|PK_Tyr_Ser-Thr|APH|HAD_\d+|Hexokinase|"
     r"Pribosyltran|NTP_transf_\d+|Guanylate_kinase|Adenylate_kin)", "EC 2.7 - phosphotransferase (kinase)"),
    # EC 2.1 — one-carbon transfer (methylation)
    (r"^(Methyltransf_\d+|MTS|Ubie_methyltran|AdoMet_MTases|FtsJ|"
     r"N6-adenineMlase|SAM_MT)", "EC 2.1 - methyltransferase"),
    # EC 2.3 / 2.4 — acyltransferases and glycosyltransferases
    (r"^(Acetyltransf_\d+|Acyl_transf_\d+|CAT_RBD|Condensation|Chal_sti_synt_C|"
     r"Glycos_transf_\d+|GT\d+|DXDXXD|DPM1_like)", "EC 2.3/2.4 - acyl- / glycosyltransferase"),
    # EC 4 — lyases
    (r"(Lyase|lyase|Synthase|synthase|Dehydratase|Dehydratase|Aldolase|Aconitase|"
     r"Fumarase|Enolase|Terpene_synth|Terpene_synth_C|ECH_\d*|Hydratase)", "EC 4 - lyase"),
    # EC 5 — isomerases
    (r"(Isomerase|isomerase|Epimerase_2|Mutase|mutase|Chalcone|"
     r"PpiC|Cyclophilin|Rotamase)", "EC 5 - isomerase"),
    # EC 6 — ligases
    (r"(Ligase|ligase|tRNA-synt_\d+|tRNA-synt_1|AMP-binding|Luciferase|"
     r"Ubiquitin-conjugating|SUMO|MoeB|ThiF)", "EC 6 - ligase"),
    # EC 7 — translocases
    (r"(MFS_\d+|Sugar_tr|Mito_carr|ABC2_membrane_\d+|Na_H_Exchanger|PQ-loop|"
     r"Transp_\w+|Permease|Major_Facilitator|Mpv17_PMP22|Gpr1_Fun34_YaaH|"
     r"Ammonium_transp|Nramp|K_trans|Cation_ATPase|E1-E2_ATPase|HATPase_c_\d)", "EC 7 - translocase"),
]
EC_ORDER = [c for _, c in EC_PATTERNS]

MF_PATTERNS = [
    (r"^Glyco_hydro", "Hydrolase activity (EC 3)"),
    (r"^Glycos_transf|^Glyco_transf|^GT[0-9]|^Mannosyl", "Transferase activity (EC 2)"),
    (r"^Polysacc_deac|^Pectinesterase|^Chitin", "Hydrolase activity (EC 3)"),
    (r"^Pec_lyase|^Lyase", "Lyase activity (EC 4)"),
    (r"^CBM_|^Lectin|^Concanavalin|^LysM", "Binding: carbohydrate"),
    (r"^Epimerase|^Aldose_epim|^Isomerase", "Isomerase activity (EC 5)"),
    (r"^Pkinase|^Kinase|^ABC1|^Choline_kin", "Transferase activity (EC 2)"),
    (r"^Peptidase|^Protease|^Subtilase|^Amidase|^Asp$", "Hydrolase activity (EC 3)"),
    (r"^tRNA-synt|^tRNA_synt|^DNA_ligase|^AMP-binding|^ATP-grasp", "Ligase activity (EC 6)"),
    (r"^ABC_tran|^ABC2_membrane|^ABC_membrane|^MFS|^Sugar_tr|^Aa_trans|^AA_permease|^Mito_carr",
     "Translocase / transporter activity (EC 7)"),
    (r"^HSP|^DnaJ|^Cpn|^Chaperone|^TCP-1", "Chaperone activity"),
    (r"^Thioredoxin|^Glutaredoxin|^GSHPx|^Peroxidase", "Antioxidant activity"),
    (r"^Actin|^Tubulin|^Myosin|^Kinesin|^Dynein", "Structural molecule activity"),
    (r"^Ank|^TPR|^WD40|^LRR|^SH3|^WW|^F-box|^UBA|^BTB", "Binding: protein"),
    (r"^RRM|^KH|^DEAD|^Helicase|^zf-CCCH", "Binding: nucleic acid"),
    (r"^HTH|^Fungal_trans|^Zn_clus|^bZIP|^zf-C2H2|^HMG|^Histone|^Bromodomain", "Transcription regulator activity"),
    (r"^EF-hand|^EF_hand|^C2$|^PH$", "Binding: small molecule or ligand"),

    # EC-level splits inside the large heterogeneous clans (NADP_Rossmann, P-loop_NTPase,
    # TIM_barrel, Beta_propeller), where the clan label alone would collapse distinct
    # enzyme classes into one.
    (r"^adh_short|^ADH_|^Aldo_ket_red|^DAO$|^Pyr_redox|^FMN_dh|^Oxidored|^IMPDH|^KR$|^FAD_binding|^NAD_binding_[0-9]+_?.*(?:dh|red)",
     "Oxidoreductase activity (EC 1)"),
    (r"^2OG-FeII_Oxy|^p450|^FMO|^GFO_IDH|^GMC_oxred|^DAO|^NmrA_", "Oxidoreductase activity (EC 1)"),
    (r"^NAD_binding|^NADP|^cofactor", "Binding: small molecule or ligand"),
    (r"^Methyltransf|^Acetyltransf|^Aminotran|^NTP_transferase|^Glu_synthase|^Pribosyltran",
     "Transferase activity (EC 2)"),
    (r"^Epimerase|^Aldose_epim", "Isomerase activity (EC 5)"),
    (r"^NmrA$", "Transcription regulator activity"),
    (r"^AAA|^Helicase|^DEAD|^ResIII|^MMR_HSR1|^Ras$|^Arf$|^Amidohydro|^SNF2", "Hydrolase activity (EC 3)"),
    (r"^NACHT|^Roc$", "Signal sensor / transducer activity"),
    (r"^SMC_N", "Binding: nucleic acid"),
    (r"^WD40|^ANAPC|^Kelch|^PD40|^Beta-prop|^Nup1", "Binding: protein"),
    (r"^eIF", "Translation regulator activity"),
    (r"^Lactonase|^DHDPS|^HpcH|^Glu_synthase", "Hydrolase activity (EC 3)"),
]

# Families that carry weight here but whose clan is uninformative and whose description
# does not hit a keyword; each is placed on the MF axis from its own definition.
MF_FAMILY_OVERRIDE = {
    "ECH_2": "Lyase activity (EC 4)", "PigL_N": "Oxidoreductase activity (EC 1)",
    "GFA": "Lyase activity (EC 4)", "HLH": "Transcription regulator activity",
    "Tyrosinase": "Oxidoreductase activity (EC 1)", "Rieske": "Oxidoreductase activity (EC 1)",
    "Ring_hydroxyl_A": "Oxidoreductase activity (EC 1)",
    "TPP_enzyme_C": "Transferase activity (EC 2)", "TPP_enzyme_N": "Transferase activity (EC 2)",
    "Biotin_lipoyl": "Transferase activity (EC 2)", "Transket_pyr": "Transferase activity (EC 2)",
    "Aldolase_II": "Lyase activity (EC 4)", "Aconitase": "Lyase activity (EC 4)",
    "zf-MYND": "Binding: protein", "FHA": "Binding: protein", "SesA": "Binding: protein",
    "VWA": "Binding: protein", "Fasciclin": "Binding: protein", "SUR7": "Binding: protein",
    "Gpr1_Fun34_YaaH": "Translocase / transporter activity (EC 7)",
    "Na_H_Exchanger": "Translocase / transporter activity (EC 7)",
    "PQ-loop": "Translocase / transporter activity (EC 7)",
    "Mpv17_PMP22": "Translocase / transporter activity (EC 7)",
    "Bac_rhamnosid_C": "Hydrolase activity (EC 3)", "Hydantoinase_A": "Hydrolase activity (EC 3)",
    "Hydant_A_N": "Hydrolase activity (EC 3)", "Arginase": "Hydrolase activity (EC 3)",
    "PAP2": "Hydrolase activity (EC 3)", "BCS1_N": "Hydrolase activity (EC 3)",
    "RGS": "Molecular function regulator", "RhoGAP": "Molecular function regulator",
    "RhoGEF": "Molecular function regulator",
    "MreB_Mbl": "Structural molecule activity",
    "G-patch": "Binding: nucleic acid", "HGTP_anticodon": "Binding: nucleic acid",
    "CBS": "Binding: small molecule or ligand", "HMA": "Binding: small molecule or ligand",
    "SPX": "Binding: small molecule or ligand", "EthD": "Catalytic activity, unspecified",
}


# ---------------------------------------------------------------- clan assignment
# Explicit classes for the clans that carry the bulk of the domains in the FSSC annotation
# (the top clans cover >90% of the gene-domain records), assigned from the clan name and
# its members' descriptions.
CLAN_MAP = {
    "NADP_Rossmann": "Oxidoreduction", "MFS": "Transport", "TPR": "Structural repeat / protein interaction",
    "Ank": "Structural repeat / protein interaction", "AB_hydrolase": "Hydrolysis",
    "TIM_barrel": "Carbohydrate metabolism and CAZymes", "Beta_propeller": "Structural repeat / protein interaction",
    "PKinase": "Signal transduction", "Fungal_trans": "Transcription and chromatin", "APC": "Transport",
    "Zn2_Cys6": "Transcription and chromatin", "Cupin": "Other enzymatic", "RING": "Protein folding, modification and turnover",
    "Acetyltrans": "Transfer of groups", "Thioredoxin": "Oxidoreduction", "HTH": "Transcription and chromatin",
    "HAD": "Hydrolysis", "PLP_aminotran": "Transfer of groups", "SH3": "Structural repeat / protein interaction",
    "GT-A": "Carbohydrate metabolism and CAZymes", "E-set": "Structural repeat / protein interaction",
    "C2H2-zf": "Transcription and chromatin", "GroES": "Protein folding, modification and turnover",
    "RRM": "RNA processing and RNA binding", "Thiolase": "Lipid and secondary metabolism",
    "F-box": "Protein folding, modification and turnover", "EF_hand": "Signal transduction",
    "ANL": "Lipid and secondary metabolism", "6_Hairpin": "Carbohydrate metabolism and CAZymes",
    "ALR_C-like": "Other enzymatic", "Peptidase_CA": "Hydrolysis", "GST_C": "Oxidoreduction",
    "ALDH-like": "Oxidoreduction", "GT-B": "Carbohydrate metabolism and CAZymes",
    "Concanavalin": "Structural repeat / protein interaction", "FAD_PCMH": "Oxidoreduction",
    "ATP-grasp": "Transfer of groups", "HUP": "Translation and ribosome biogenesis",
    "OB": "Translation and ribosome biogenesis", "Peptidase_MA": "Hydrolysis", "ABC-2": "Transport",
    "Pec_lyase-like": "Carbohydrate metabolism and CAZymes", "PH": "Signal transduction",
    "NAD_Ferredoxin": "Oxidoreduction", "Ubiquitin": "Protein folding, modification and turnover",
    "Alk_phosphatase": "Hydrolysis", "CoA-acyltrans": "Transfer of groups",
    "GBD": "Structural repeat / protein interaction", "ZnExoPePases": "Hydrolysis", "Peptidase_AA": "Hydrolysis",
    "Calcineurin": "Hydrolysis", "NTN": "Protein folding, modification and turnover",
    "LTE-like": "Transport", "His_Kinase_A": "Signal transduction", "DMT": "Transport",
    "CU_oxidase": "Oxidoreduction", "Glutaminase_I": "Transfer of groups", "HotDog": "Lipid and secondary metabolism",
    "PAS_Fold": "Signal transduction", "6PGD_C": "Oxidoreduction", "Histone": "Transcription and chromatin",
    "ABC_membrane": "Transport", "Zn_Beta_Ribbon": "Structural repeat / protein interaction",
    "Phosphatase": "Hydrolysis", "Patatin": "Lipid and secondary metabolism", "Metallo-HOrase": "Hydrolysis",
    "NTF2": "Trafficking and secretion", "LRR": "Structural repeat / protein interaction",
    "FAD-oxidase_C": "Oxidoreduction", "B_Fructosidase": "Carbohydrate metabolism and CAZymes",
    "SGNH_hydrolase": "Hydrolysis", "AMP-binding_C": "Lipid and secondary metabolism",
    "GPCR_A": "Signal transduction", "PP-binding": "Lipid and secondary metabolism",
    "Gal_mutarotase": "Carbohydrate metabolism and CAZymes", "Leu-IlvD": "Oxidoreduction",
    "Beta-lactamase": "Hydrolysis", "bZIP": "Transcription and chromatin",
    "tRNA_synt_II": "Translation and ribosome biogenesis", "UBC": "Protein folding, modification and turnover",
    "Acyl-CoA_dh": "Lipid and secondary metabolism", "Cyclin": "Transcription and chromatin",
    "BRCT-like": "DNA replication, repair and recombination", "Chaperone-J": "Protein folding, modification and turnover",
    "RNase_H": "DNA replication, repair and recombination", "2heme_cytochrom": "Oxidoreduction",
    "His_phosphatase": "Hydrolysis", "POZ": "Protein folding, modification and turnover",
    "UBA": "Protein folding, modification and turnover", "CheY": "Signal transduction",
    "4Fe-4S": "Oxidoreduction", "FMN-binding": "Oxidoreduction", "Peptidase_CD": "Hydrolysis",
    "GADPH_aa-bio_dh": "Oxidoreduction", "Steroid_dh": "Lipid and secondary metabolism",
    "PAN": "Structural repeat / protein interaction", "SAF": "Transcription and chromatin",
    "PRTase-like": "Transfer of groups", "Peptidase_ME": "Hydrolysis", "PF": "Cytoskeleton and motility",
    "NUDIX": "Hydrolysis", "Golgi-transport": "Trafficking and secretion", "SNARE-fusion": "Trafficking and secretion",
    "PBP": "Transport", "Flavoprotein": "Oxidoreduction", "PA14": "Structural repeat / protein interaction",
    "CDA": "DNA replication, repair and recombination", "Arrestin_N-like": "Trafficking and secretion",
    "C2": "Signal transduction", "GHD": "Carbohydrate metabolism and CAZymes", "DHQS": "Oxidoreduction",
    "GT-C": "Cell wall and membrane biogenesis", "Transthyretin": "Other enzymatic", "KH": "RNA processing and RNA binding",
    "PDDEXK": "Other enzymatic", "NTP_transf": "Transfer of groups", "FUSC": "Transport",
    "Acyl-CoA_desat": "Lipid and secondary metabolism",
    "P-loop_NTPase": "Transfer of groups",   # mostly ATP/GTP-binding cassettes; family overrides below
}

# Families whose job differs from the rest of their clan.
FAMILY_OVERRIDE = {
    "ABC_tran": "Transport", "ABC2_membrane": "Transport", "ABC_membrane": "Transport",
    "AAA": "Protein folding, modification and turnover",
    "Sugar_tr": "Transport", "MFS_1": "Transport", "MFS_4": "Transport",
    "ADH_N": "Oxidoreduction", "ADH_N_2": "Oxidoreduction", "ADH_short": "Oxidoreduction",
    "ADH_zinc_N": "Oxidoreduction", "Aldedh": "Oxidoreduction", "Aldo_ket_red": "Oxidoreduction",
    "DAO": "Oxidoreduction", "FMO-like": "Oxidoreduction", "GFO_IDH_MocA": "Oxidoreduction",
    "FAD_binding_2": "Oxidoreduction", "FAD_binding_3": "Oxidoreduction", "FAD_binding_4": "Oxidoreduction",
    "NAD_binding_8": "Oxidoreduction", "NAD_binding_10": "Oxidoreduction",
    "Ribosomal_L3": "Translation and ribosome biogenesis", "Ribosomal_L16": "Translation and ribosome biogenesis",
    "HSP70": "Protein folding, modification and turnover", "Hexokinase_1": "Carbohydrate metabolism and CAZymes",
    "ECH_1": "Lipid and secondary metabolism", "Cellulase": "Carbohydrate metabolism and CAZymes",
    "Glyco_hydro_47": "Carbohydrate metabolism and CAZymes", "Glyco_hydro_43": "Carbohydrate metabolism and CAZymes",
    "Glyco_hydro_28": "Carbohydrate metabolism and CAZymes", "Glyco_hydro_18": "Carbohydrate metabolism and CAZymes",
    "Epimerase": "Carbohydrate metabolism and CAZymes", "NmrA": "Transcription and chromatin",
    "Chitin_synth_2": "Cell wall and membrane biogenesis", "Chitin_synth_1": "Cell wall and membrane biogenesis",
    "DEAD": "RNA processing and RNA binding", "SNF2_N": "Transcription and chromatin",
    "RTA1": "Transport", "Mito_carr": "Transport", "Aa_trans": "Transport", "AA_permease": "Transport",
    "Fungal_trans": "Transcription and chromatin", "Zn_clus": "Transcription and chromatin",
    "AMP-binding": "Lipid and secondary metabolism", "NACHT": "Signal transduction",
    "WD40": "Structural repeat / protein interaction", "Ank_2": "Structural repeat / protein interaction",
    "HET": "Other enzymatic",  # heterokaryon incompatibility, fungal-specific
}


# Families WITHOUT a clan (13,524 of the 25,545 in Pfam 38.0, carrying 14% of this project's
# domain records) cannot inherit a clan class, so each is assigned from its own name and
# description. The entries below are the no-clan families that actually carry weight here.
NOCLAN_OVERRIDE = {
    "HET": "Other enzymatic", "p450": "Oxidoreduction", "Mito_carr": "Transport",
    "Amidase": "Hydrolysis", "E1-E2_ATPase": "Transport", "Peptidase_S8": "Hydrolysis",
    "Cyt-b5": "Oxidoreduction", "Glyco_hydro_3_C": "Carbohydrate metabolism and CAZymes",
    "Cation_ATPase_N": "Transport", "Cation_ATPase_C": "Transport", "WSC": "Cell wall and membrane biogenesis",
    "OPT": "Transport", "CN_hydrolase": "Hydrolysis", "HeLo": "Unknown function",
    "Helo_like_N": "Lipid and secondary metabolism", "CFEM": "Structural repeat / protein interaction",
    "M20_dimer": "Hydrolysis", "ABC_trans_N": "Transport", "Clr5": "Unknown function",
    "CBM_1": "Carbohydrate metabolism and CAZymes", "Chitin_bind_1": "Carbohydrate metabolism and CAZymes",
    "LysM": "Carbohydrate metabolism and CAZymes", "Polysacc_deac_1": "Carbohydrate metabolism and CAZymes",
    "Glyco_hydro_3": "Carbohydrate metabolism and CAZymes", "Glyco_hydro_5": "Carbohydrate metabolism and CAZymes",
    "Glyco_hydro_16": "Carbohydrate metabolism and CAZymes", "Glyco_hydro_17": "Carbohydrate metabolism and CAZymes",
    "Cpn60_TCP1": "Protein folding, modification and turnover", "Bromodomain": "Transcription and chromatin",
    "HMG_box": "Transcription and chromatin", "HLH": "Transcription and chromatin", "zf-MYND": "Transcription and chromatin",
    "zf-CCCH": "RNA processing and RNA binding", "PX": "Signal transduction", "WW": "Structural repeat / protein interaction",
    "FHA": "Signal transduction", "RGS": "Signal transduction", "RhoGAP": "Signal transduction",
    "Tyrosinase": "Oxidoreduction", "Rieske": "Oxidoreduction", "Ring_hydroxyl_A": "Oxidoreduction",
    "GFA": "Oxidoreduction", "ECH_2": "Lipid and secondary metabolism", "OrsD": "Lipid and secondary metabolism",
    "Isochorismatase": "Hydrolysis", "Hydantoinase_A": "Hydrolysis", "Hydant_A_N": "Hydrolysis",
    "TPP_enzyme_C": "Transfer of groups", "TPP_enzyme_N": "Transfer of groups", "Biotin_lipoyl": "Transfer of groups",
    "Transket_pyr": "Carbohydrate metabolism and CAZymes", "Aldolase_II": "Carbohydrate metabolism and CAZymes",
    "Na_Ca_ex": "Transport", "Na_H_Exchanger": "Transport", "Cation_efflux": "Transport",
    "PQ-loop": "Transport", "Grp1_Fun34_YaaH": "Transport", "Dynamin_M": "Trafficking and secretion",
    "MreB_Mbl": "Cytoskeleton and motility", "Bac_rhamnosid_C": "Carbohydrate metabolism and CAZymes",
    "SesA": "Other enzymatic", "OrsD_x": "Lipid and secondary metabolism",
}

# Fallback: match the description when the clan is not in CLAN_MAP.
KEYWORD_RULES = [
    ("Carbohydrate metabolism and CAZymes", ["glycosyl hydrolase", "glycosyltransferase", "glycoside hydrolase",
                                             "carbohydrate", "glucan synthase", "chitin synth", "sugar",
                                             "epimerase", "cellulase", "amylase", "mannosidase"]),
    ("Transport", ["transport", "permease", "symport", "antiport", "channel", "carrier"]),
    ("Transcription and chromatin", ["transcription factor", "transcriptional", "histone", "homeobox", "chromatin"]),
    ("Translation and ribosome biogenesis", ["ribosomal", "trna synthetase", "translation", "elongation factor"]),
    ("RNA processing and RNA binding", ["rna-binding", "rna recognition", "rna helicase", "splic", "rnase"]),
    ("DNA replication, repair and recombination", ["dna polymerase", "dna repair", "recombin", "topoisomerase",
                                                   "transposase", "integrase", "dna ligase"]),
    ("Protein folding, modification and turnover", ["chaperone", "heat shock", "ubiquitin", "proteasome",
                                                    "peptidyl-prolyl", "disulfide"]),
    ("Cytoskeleton and motility", ["actin", "tubulin", "myosin", "kinesin", "dynein", "cytoskelet", "profilin"]),
    ("Signal transduction", ["kinase", "gtpase", "g-protein", "calmodulin", "ef hand", "second messenger"]),
    ("Oxidoreduction", ["dehydrogenase", "oxidase", "oxygenase", "reductase", "peroxidase", "catalase",
                        "cytochrome", "oxidoreductase", "flavoprotein", "ferredoxin"]),
    ("Hydrolysis", ["hydrolase", "protease", "peptidase", "esterase", "lipase", "phosphatase", "nuclease",
                    "amidase", "lactamase", "esterase"]),
    ("Transfer of groups", ["transferase", "synthase", "phosphorylase", "ligase"]),
    ("Lipid and secondary metabolism", ["lipid", "fatty acid", "sterol", "polyketide", "terpene", "alkaloid",
                                        "phospholipid"]),
    ("Trafficking and secretion", ["snare", "coatomer", "rab ", "vacuolar sorting", "nuclear transport"]),
    ("Structural repeat / protein interaction", ["ankyrin", "tetratricopeptide", "wd domain", "leucine rich repeat",
                                                 "immunoglobulin", "lectin", "armadillo"]),
]
UNKNOWN_MARKERS = ("unknown function", "uncharacteri", "duf", "putative uncharacterized", "protein of unknown")


def describe(acc):
    """Full description from pfamA.txt (column 4)."""
    with gzip.open(PfamA, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] == acc:
                return p[3]
    return ""


def main():
    rows = []
    with gzip.open(CLANS, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            acc, clan_id, clan_name, name, desc = p[0], p[1], p[2], p[3], p[4]
            if clan_name.strip() == "":      # no clan in Pfam 38.0 -> must go by definition
                clan_id, clan_name = "", ""
            low = desc.lower()
            clan_id, clan_name, name, desc = clan_id.strip(), clan_name.strip(), name.strip(), desc.strip()
            if name in FAMILY_OVERRIDE:
                klass, method = FAMILY_OVERRIDE[name], "clan_override"
            elif name in NOCLAN_OVERRIDE and clan_name == "":
                klass, method = NOCLAN_OVERRIDE[name], "definition"
            elif clan_name in CLAN_MAP:
                klass, method = CLAN_MAP[clan_name], "clan"
            elif any(m in low for m in UNKNOWN_MARKERS):
                klass, method = "Unknown function", "unknown"
            else:
                klass, method = None, None
                hay = (desc + " " + name).lower()
                for k, kws in KEYWORD_RULES:
                    if any(w in hay for w in kws):
                        klass, method = k, "keyword"
                        break
                if klass is None:
                    klass, method = "No assignment", "unmatched"


            mf = MF_FAMILY_OVERRIDE.get(name)
            for pat, mcls in MF_PATTERNS if mf is None else []:
                if _re.match(pat, name):
                    mf = mcls
                    break
            if mf is None:
                mf = SEMANTIC_TO_MF.get(klass, "No assignment")
            ec = "non-enzymatic or unassigned"
            if _re.search(r"\(EC \d\)", mf):
                for pat, ecls in EC_PATTERNS:
                    if _re.match(pat, name) or _re.search(pat, desc):
                        ec = ecls
                        break
                else:
                    ec = "EC %s - other at this level" % _re.search(r"\(EC (\d)\)", mf).group(1)
            rows.append(dict(pfam_acc=acc, pfam_name=name, description=desc, clan_id=clan_id,
                             clan_name=clan_name, functional_class=klass,
                             go_mf_class=mf, ec_subclass=ec, method=method))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(OUT_CLASSES, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["axis", "class", "definition"])
        w.writerows([("pathway-style (internal)", k, v) for k, v in CLASSES])
        w.writerows([("GO MF / EC (manuscript)", k, v) for k, v in MF_CLASSES])
    total = len(rows)
    print(f"families classified: {total:,}")
    from collections import Counter
    print("\nmethod:", dict(Counter(r["method"] for r in rows)))
    print("\nGO MF axis distribution (% of families):")
    cnt = Counter(r["go_mf_class"] for r in rows)
    for k, v in cnt.most_common():
        print(f"  {k:44s} {v:6,d}  {v / total * 100:5.1f}%")
    print(f"\nwrote {os.path.relpath(OUT, ROOT)} and {os.path.relpath(OUT_CLASSES, ROOT)}")


if __name__ == "__main__":
    main()
