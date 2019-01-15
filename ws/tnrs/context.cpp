#include "context.h"

using std::string;
using std::map;

AllContexts::AllContexts(const std::initializer_list<Context>& contexts)
{
    for(auto& context: contexts)
        insert({context.name,context});
}


// NOTE: These context definitions were taken from `taxomachine/src/main/java/org/opentree/taxonomy/contexts/ContextDescription.java`
//       on Jan 15, 2019.

AllContexts all_contexts = {
    // Name             Group                         Index suffix        Node name string    ott id      Code
    {"All life",        /* ContextGroup.LIFE, */      "",                 "life",             805080L,    Nomenclature::Undefined},

    // MICROBES group
    {"Bacteria",        /* ContextGroup.MICROBES, */  "Bacteria",         "Bacteria",         844192L,    Nomenclature::ICNP},
    {"SAR group",       /* ContextGroup.MICROBES, */  "SAR",              "SAR",              5246039L,   Nomenclature::Undefined},
    {"Archaea",         /* ContextGroup.MICROBES, */  "Archaea",          "Archaea",          996421L,    Nomenclature::ICNP},
    {"Excavata",        /* ContextGroup.MICROBES, */  "Excavata",         "Excavata",         2927065L,   Nomenclature::Undefined},
    {"Amoebozoa",       /* ContextGroup.MICROBES, */  "Amoebae",          "Amoebozoa",        1064655L,   Nomenclature::ICZN},
    {"Centrohelida",    /* ContextGroup.MICROBES, */  "Centrohelida",     "Centrohelida",     755852L,    Nomenclature::ICZN},
    {"Haptophyta",      /* ContextGroup.MICROBES, */  "Haptophyta",       "Haptophyta",       151014L,    Nomenclature::Undefined},
    {"Apusozoa",        /* ContextGroup.MICROBES, */  "Apusozoa",         "Apusozoa",         671092L,    Nomenclature::ICZN},
    {"Diatoms",         /* ContextGroup.MICROBES, */  "Diatoms",          "Bacillariophyta",  5342311L,   Nomenclature::ICN},
    {"Ciliates",        /* ContextGroup.MICROBES, */  "Ciliates",         "Ciliophora",       302424L,    Nomenclature::Undefined},
    {"Forams",          /* ContextGroup.MICROBES, */  "Forams",           "Foraminifera",     936399L,    Nomenclature::ICZN},

    // ANIMALS group
    {"Animals",         /* ContextGroup.ANIMALS, */   "Animals",          "Metazoa",          691846L,    Nomenclature::ICZN},
    {"Birds",           /* ContextGroup.ANIMALS, */   "Birds",            "Aves",             81461L,     Nomenclature::ICZN},
    {"Tetrapods",       /* ContextGroup.ANIMALS, */   "Tetrapods",        "Tetrapoda",        229562L,    Nomenclature::ICZN},
    {"Mammals",         /* ContextGroup.ANIMALS, */   "Mammals",          "Mammalia",         244265L,    Nomenclature::ICZN},
    {"Amphibians",      /* ContextGroup.ANIMALS, */   "Amphibians",       "Amphibia",         544595L,    Nomenclature::ICZN},
    {"Vertebrates",     /* ContextGroup.ANIMALS, */   "Vertebrates",      "Vertebrata",       801601L,    Nomenclature::ICZN},
    {"Arthropods",      /* ContextGroup.ANIMALS, */   "Arthopods",        "Arthropoda",       632179L,    Nomenclature::ICZN},
    {"Molluscs",        /* ContextGroup.ANIMALS, */   "Molluscs",         "Mollusca",         802117L,    Nomenclature::ICZN},
    {"Nematodes",       /* ContextGroup.ANIMALS, */   "Nematodes",        "Nematoda",         395057L,    Nomenclature::ICZN},
    {"Platyhelminthes", /* ContextGroup.ANIMALS, */   "Platyhelminthes",  "Platyhelminthes",  555379L,    Nomenclature::ICZN},
    {"Annelids",        /* ContextGroup.ANIMALS, */   "Annelids",         "Annelida",         941620L,    Nomenclature::ICZN},
    {"Cnidarians",      /* ContextGroup.ANIMALS, */   "Cnidarians",       "Cnidaria",         641033L,    Nomenclature::ICZN},
    {"Arachnids",       /* ContextGroup.ANIMALS, */   "Arachnids",        "Arachnida",        511967L,    Nomenclature::ICZN},
    {"Insects",         /* ContextGroup.ANIMALS, */   "Insects",          "Insecta",          1062253L,   Nomenclature::ICZN},

    // FUNGI group
    {"Fungi",           /* ContextGroup.FUNGI, */     "Fungi",            "Fungi",            352914L,    Nomenclature::ICN},
    {"Basidiomycetes",  /* ContextGroup.FUNGI, */     "Basidiomycetes",   "Basidiomycota",    634628L,    Nomenclature::ICN},
    {"Ascomycetes",     /* ContextGroup.FUNGI, */     "Ascomycota",       "Ascomycota",       439373L,    Nomenclature::ICN},
    
    // PLANTS group
    {"Land plants",     /* ContextGroup.PLANTS, */    "Plants",           "Embryophyta",      5342313L,   Nomenclature::ICN},
    {"Hornworts",       /* ContextGroup.PLANTS, */    "Anthocerotophyta", "Anthocerotophyta", 738980L,    Nomenclature::ICN},
    {"Mosses",          /* ContextGroup.PLANTS, */    "Bryophyta",        "Bryophyta",        246594L,    Nomenclature::ICN},
    {"Liverworts",      /* ContextGroup.PLANTS, */    "Marchantiophyta",  "Marchantiophyta",  56601L,     Nomenclature::ICN},
    {"Vascular plants", /* ContextGroup.PLANTS, */    "Tracheophyta",     "Tracheophyta",     10210L,     Nomenclature::ICN},
    {"Club mosses",     /* ContextGroup.PLANTS, */    "Lycopodiophyta",   "Lycopodiophyta",   144803L,    Nomenclature::ICN},
    {"Ferns",           /* ContextGroup.PLANTS, */    "Moniliformopses",  "Moniliformopses",  166292L,    Nomenclature::ICN},
    {"Seed plants",     /* ContextGroup.PLANTS, */    "Spermatophyta",    "Spermatophyta",    10218L,     Nomenclature::ICN},
    {"Flowering plants",/* ContextGroup.PLANTS, */    "Magnoliophyta",    "Magnoliophyta",    99252L,     Nomenclature::ICN},
    {"Monocots",        /* ContextGroup.PLANTS, */    "Monocots",         "Liliopsida",       1058517L,   Nomenclature::ICN},
    {"Eudicots",        /* ContextGroup.PLANTS, */    "Eudicots",         "eudicotyledons",   431495L,    Nomenclature::ICN},
    {"Rosids",          /* ContextGroup.PLANTS, */    "Rosids",           "rosids",           1008296L,   Nomenclature::ICN},
    {"Asterids",        /* ContextGroup.PLANTS, */    "Asterids",         "asterids",         1008294L,   Nomenclature::ICN},
    {"Asterales",       /* ContextGroup.PLANTS, */    "Asterales",        "Asterales",        1042120L,   Nomenclature::ICN},
    {"Asteraceae",      /* ContextGroup.PLANTS, */    "Asteraceae",       "Asteraceae",       46248L,     Nomenclature::ICN},    
    {"Aster",           /* ContextGroup.PLANTS, */    "Aster",            "Aster",            409712L,    Nomenclature::ICN},    
    {"Symphyotrichum",  /* ContextGroup.PLANTS, */    "Symphyotrichum",   "Symphyotrichum",   1058735L,   Nomenclature::ICN},
    {"Campanulaceae",   /* ContextGroup.PLANTS, */    "Campanulaceae",    "Campanulaceae",    1086303L,   Nomenclature::ICN},
    {"Lobelia",         /* ContextGroup.PLANTS, */    "Lobelia",          "Lobelia",          1086294L,   Nomenclature::ICN}
};
