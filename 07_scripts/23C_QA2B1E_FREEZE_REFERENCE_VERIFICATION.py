#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib

ROOT = Path.cwd()

DOCX = (
    ROOT /
    "04_results/manuscript_update/Step23/"
    "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA2B1C_QA2E_LOCKED_EDITS.docx"
)

META = ROOT / "00_metadata/Step23_manuscript_update"
RESULT = ROOT / "04_results/manuscript_update/Step23"

OUT_CSV = META / "STEP23C_QA2B1E_REFERENCE_VERIFICATION_LOCK.csv"
OUT_TXT = RESULT / "STEP23C_QA2B1E_REFERENCE_VERIFICATION_DOSSIER.txt"
OUT_SHA = META / "STEP23C_QA2B1E_REFERENCE_VERIFICATION_CHECKSUMS.sha256"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


if not DOCX.exists():
    raise SystemExit(f"STOP: missing QA2B1C DOCX:\n{DOCX}")

docx_sha_before = sha256(DOCX)

refs = [
    {
        "accession": "GSE197677",
        "first_author": "Sho Okuda",
        "authors": (
            "Sho Okuda; Kenoki Ohuchida; Shoichi Nakamura; "
            "Chikanori Tsutsumi; Kyoko Hisano; Yuki Mochida; "
            "Jun Kawata; Yoshiki Ohtsubo; Tomohiko Shinkawa; "
            "Chika Iwamoto; Nobuhiro Torata; Yusuke Mizuuchi; "
            "Koji Shindo; Taiki Moriyama; Kohei Nakata; "
            "Takehiro Torisu; Takashi Morisaki; Takanari Kitazono; "
            "Yoshinao Oda; Masafumi Nakamura"
        ),
        "title": (
            "Neoadjuvant chemotherapy enhances anti-tumor immune response "
            "of tumor microenvironment in human esophageal squamous cell carcinoma"
        ),
        "journal": "iScience",
        "year": "2023",
        "volume": "26",
        "issue": "4",
        "article": "106480",
        "doi": "10.1016/j.isci.2023.106480",
        "pmid": "37091252",
        "dataset_relation": "GEO GSE197677 / BioProject PRJNA811547",
        "verification_status": "VERIFIED_PRIMARY_PUBLICATION",
    },
    {
        "accession": "GSE221561",
        "first_author": "Yushang Yang",
        "authors": (
            "Yushang Yang; Yanguo Li; Haopeng Yu; Zhenyu Ding; "
            "Longqi Chen; Xiaoxi Zeng; Shunmin He; Qi Liao; "
            "Yi Zhao; Yong Yuan"
        ),
        "title": (
            "Comprehensive landscape of resistance mechanisms for "
            "neoadjuvant therapy in esophageal squamous cell carcinoma "
            "by single-cell transcriptomics"
        ),
        "journal": "Signal Transduction and Targeted Therapy",
        "year": "2023",
        "volume": "8",
        "issue": "1",
        "article": "298",
        "doi": "10.1038/s41392-023-01518-0",
        "pmid": "37563120",
        "dataset_relation": "GEO GSE221561",
        "verification_status": "VERIFIED_PRIMARY_PUBLICATION",
    },
    {
        "accession": "OMIX005710",
        "first_author": "Gang Ji",
        "authors": (
            "Gang Ji; Qi Yang; Song Wang; Xiaolong Yan; Qiuxiang Ou; "
            "Li Gong; Jinbo Zhao; Yongan Zhou; Feng Tian; Jie Lei; "
            "Xiaorong Mu; Jian Wang; Tao Wang; Xiaoping Wang; "
            "Jianyong Sun; Jipeng Zhang; Chenghui Jia; Tao Jiang; "
            "Ming-Gao Zhao; Qiang Lu"
        ),
        "title": (
            "Single-cell profiling of response to neoadjuvant "
            "chemo-immunotherapy in surgically resectable esophageal "
            "squamous cell carcinoma"
        ),
        "journal": "Genome Medicine",
        "year": "2024",
        "volume": "16",
        "issue": "1",
        "article": "49",
        "doi": "10.1186/s13073-024-01320-9",
        "pmid": "38566201",
        "dataset_relation": "OMIX005710 / BioProject PRJCA016745",
        "verification_status": "VERIFIED_PRIMARY_PUBLICATION",
    },
]

META.mkdir(parents=True, exist_ok=True)
RESULT.mkdir(parents=True, exist_ok=True)

fields = list(refs[0].keys())

with open(
    OUT_CSV,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(refs)


lines = [
    "STEP23C-QA2B1E REFERENCE VERIFICATION DOSSIER",
    "=" * 80,
    "",
    "STATUS: READ_ONLY_BIBLIOGRAPHIC_VERIFICATION",
    "",
    f"QA2B1C DOCX SHA256: {docx_sha_before}",
    "",
]

for n, r in enumerate(refs, 1):
    lines += [
        "=" * 80,
        f"{n}. {r['accession']}",
        "=" * 80,
        "",
        f"Authors: {r['authors']}",
        "",
        f"Title: {r['title']}",
        f"Journal: {r['journal']}",
        f"Year: {r['year']}",
        f"Volume: {r['volume']}",
        f"Issue: {r['issue']}",
        f"Article: {r['article']}",
        f"DOI: {r['doi']}",
        f"PMID: {r['pmid']}",
        f"Dataset relation: {r['dataset_relation']}",
        f"Verification: {r['verification_status']}",
        "",
    ]

lines += [
    "=" * 80,
    "DECISION",
    "=" * 80,
    "",
    "QA2C-054 bibliographic source requirement:",
    "VERIFIED FOR GSE197677 / GSE221561 / OMIX005710",
    "",
    "IMPORTANT:",
    "This step does NOT claim that the manuscript's complete literature",
    "bibliography has been assembled. It resolves only the frozen QA2C-054",
    "source-publication verification hold.",
    "",
    "DOCX MODIFIED: NO",
    "",
    "STATUS:",
    "STEP23C_QA2B1E_REFERENCE_VERIFICATION_LOCK_COMPLETE",
    "",
]

OUT_TXT.write_text(
    "\n".join(lines),
    encoding="utf-8"
)

with open(OUT_SHA, "w", encoding="utf-8") as f:
    for p in [DOCX, OUT_CSV, OUT_TXT]:
        f.write(
            f"{sha256(p)}  {p.relative_to(ROOT)}\n"
        )

if sha256(DOCX) != docx_sha_before:
    raise SystemExit(
        "CRITICAL STOP: QA2B1C DOCX changed during verification."
    )

print()
print("=" * 88)
print("STEP23C-QA2B1E COMPLETE")
print("=" * 88)
print()
print("QA2B1C DOCX MODIFIED        : NO")
print("REFERENCES VERIFIED         : 3 / 3")
print("GSE197677                   : VERIFIED")
print("GSE221561                   : VERIFIED")
print("OMIX005710                  : VERIFIED")
print()
print("OUTPUTS:")
print(" ", OUT_CSV.relative_to(ROOT))
print(" ", OUT_TXT.relative_to(ROOT))
print(" ", OUT_SHA.relative_to(ROOT))
print()
print("STATUS:")
print("STEP23C_QA2B1E_REFERENCE_VERIFICATION_LOCK_COMPLETE")
print()
print("NEXT:")
print(
    "INSERT THE THREE VERIFIED SOURCE PUBLICATIONS INTO A NEW DOCX COPY; "
    "DO NOT MARK THE ENTIRE LITERATURE BIBLIOGRAPHY COMPLETE."
)
print()

