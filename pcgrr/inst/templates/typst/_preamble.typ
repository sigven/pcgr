// =============================================================================
// PCGR Typst preamble
// Defines brand colors, page chrome, show rules, and reusable components.
// Included via `include-before-body` in pcgr_report_pdf.qmd. Layout and
// components mirror the CPSR PDF report; colors follow the PCGR report theme
// (purple: tumor-control, blue: tumor-only - see typst/_theme.typ).
// =============================================================================

#import "typst/_theme.typ": pcgr-main

// ---------------------------------------------------------------------------
// Brand colors - derived from the main report color
// ---------------------------------------------------------------------------
#let pcgr-main-dark = pcgr-main.darken(30%)  // table headers, cover headings
#let pcgr-dark    = rgb("#2c313c")  // headings, body chrome (color_palette$bg_dark)
#let pcgr-neutral = rgb("#868686")  // empty / no-finding state (color_palette$none)
#let pcgr-rule    = color.mix((pcgr-main, 18%), (luma(205), 82%))  // hairlines
#let pcgr-stripe  = pcgr-main.lighten(94%)  // zebra rows
#let pcgr-muted   = rgb("#5f5966")  // secondary text (AA contrast on white)

// Page margins - must match the `margin` block in pcgr_report_pdf.qmd; used
// to let the cover band run full-bleed
#let pcgr-mx = 1.5cm
#let pcgr-bar = 0.32cm  // height of the full-bleed top bar (cover)

// Logo (copied together with the templates to the render directory)
#let pcgr-logo(height) = image("typst/pcgr_logo.png", height: height)

// Report-level values set by the cover, read by header/footer on later pages
#let pcgr-sample = state("pcgr-sample", "")
#let pcgr-ver    = state("pcgr-ver", "")
#let pcgr-today  = datetime.today().display("[day] [month repr:long] [year]")

// Actionability tiers (mirrors color_palette$tier_sensitivity /
// tier_resistance in R)
#let pcgr-col-tier1     = rgb("#1565C0")
#let pcgr-col-tier2     = rgb("#42A5F5")
#let pcgr-col-tier3     = rgb("#868686")
#let pcgr-col-res-tier1 = rgb("#C62828")
#let pcgr-col-res-tier2 = rgb("#E57373")

// Oncogenicity scale (mirrors color_palette$oncogenicity in R)
#let pcgr-col-onc     = rgb("#9E0142")
#let pcgr-col-lonc    = rgb("#D53E4F")
#let pcgr-col-vus     = rgb("#2c313c")
#let pcgr-col-lbenign = rgb("#78C679")
#let pcgr-col-benign  = rgb("#077009")

// ---------------------------------------------------------------------------
// Page chrome - full-bleed bands (main report color) frame every page. Top
// band (page 2+) names the section(s) on the page and carries sample ID +
// date; bottom band closes the page; footer (white) holds logo,
// disclaimer/links and page count. Margins/paper come from Quarto YAML
// ---------------------------------------------------------------------------
#let pcgr-band-top = 0.95cm   // top band, pages 2+
#let pcgr-band-bot = 0.55cm   // bottom band, all pages

// Section label for the current page: the level-1 heading(s) on the page; if
// the first heading starts below the top part of the page, the section
// continued from the previous page is listed first. Call inside `context`.
#let pcgr-section-label() = {
  let pg = here().page()
  let hs = query(heading.where(level: 1))
  let on = hs.filter(h => h.location().page() == pg)
  let before = hs.filter(h => h.location().page() < pg)
  let names = ()
  if before.len() > 0 and (
      on.len() == 0 or on.first().location().position().y > 5cm) {
    names.push(before.last().body)
  }
  for h in on { names.push(h.body) }
  if names.len() == 0 { none } else { names.join([ · ]) }
}

#set page(
  background: context {
    let pg = here().page()
    if pg == 1 {
      // cover: thin bar, the cover has its own large logo and band
      place(top, rect(width: 100%, height: pcgr-bar, fill: pcgr-main, stroke: none))
    } else {
      place(top, rect(
        width: 100%, height: pcgr-band-top, fill: pcgr-main, stroke: none,
        inset: (x: pcgr-mx),
      )[
        #set text(fill: white, size: 8.5pt)
        #align(horizon, grid(
          columns: (1fr, auto),
          align: (left + horizon, right + horizon),
          text(weight: "bold", size: 9pt, tracking: 0.1em)[#upper(pcgr-section-label())],
          [#text(weight: "bold")[#pcgr-sample.final()] · #pcgr-today],
        ))
      ])
    }
    place(bottom, rect(
      width: 100%, height: pcgr-band-bot, fill: pcgr-main, stroke: none))
  },
  header: none,
  footer: context {
    line(length: 100%, stroke: 0.5pt + pcgr-rule)
    v(3pt)
    set text(size: 7.5pt, fill: pcgr-muted)
    grid(
      columns: (auto, 1fr, auto),
      column-gutter: 8pt,
      align: (left + horizon, left + horizon, right + horizon),
      // cover page already carries the large logo
      if counter(page).get().first() == 1 { none } else { pcgr-logo(0.8cm) },
      [For research use only · pcgr v#pcgr-ver.final() ·
        #link("https://sigven.github.io/pcgr")[sigven.github.io/pcgr]],
      counter(page).display("1 of 1", both: true),
    )
  },
  // footer kept clear of the bottom band
  footer-descent: 8%,
)

// ---------------------------------------------------------------------------
// Typography
// ---------------------------------------------------------------------------
// Sans-serif font in line with the HTML report - Source Sans Pro is shipped
// with the rmarkdown R package, DejaVu Sans with the conda environment
// (fallback for symbols not covered by Source Sans Pro); both are made
// available to Typst through TYPST_FONT_PATHS (see pcgrr::write_report_pdf)
#set text(
  font: ("Source Sans Pro", "DejaVu Sans"),
  size: 10pt, fill: pcgr-dark)
#set par(leading: 0.65em)
// tables (wrapped in figures by Quarto) may break across pages
#show figure: set block(breakable: true)
// figure captions (legends): set apart from the body text - smaller, muted,
// left-aligned, with the figure number in bold
#show figure.caption: it => block(
  width: 100%, above: 0.6em, below: 1.2em,
  inset: (left: 7pt, y: 2pt),
  stroke: (left: 1.5pt + pcgr-rule),
)[
  #set align(left)
  #set par(justify: false, leading: 0.5em)
  #set text(size: 8pt, fill: pcgr-muted)
  #text(weight: "bold", fill: pcgr-dark)[
    #it.supplement #context it.counter.display(it.numbering)#it.separator]#it.body
]
// table cells: ragged-right text without hyphenation (narrow columns)
#show table.cell: set par(justify: false)
#show table.cell: set text(hyphenate: false)
// tables set slightly smaller than body text (10pt)
#show table: set text(size: 8.5pt)
// links in the main report color (clickable in the PDF)
#show link: set text(fill: pcgr-main)

// ---------------------------------------------------------------------------
// Heading show rules
// ---------------------------------------------------------------------------
#show heading.where(level: 1): it => block(above: 1.6em, below: 0.6em)[
  #set text(size: 13pt, weight: "bold", fill: pcgr-main)
  #it.body
  #v(-4pt)
  #line(length: 100%, stroke: 1.2pt + pcgr-main)
]

#show heading.where(level: 2): it => block(above: 1.2em, below: 0.4em)[
  #set text(size: 11pt, weight: "bold", fill: pcgr-dark)
  #it.body
]

#show heading.where(level: 3): it => block(above: 1em, below: 0.3em)[
  #set text(size: 10pt, weight: "bold", fill: pcgr-dark)
  #it.body
]

// ---------------------------------------------------------------------------
// Reusable components
// ---------------------------------------------------------------------------

// Tier / oncogenicity / evidence badge (inline pill)
#let pcgr-badge(label, color: pcgr-neutral) = box(
  fill: color,
  inset: (x: 5pt, y: 2pt),
  radius: 3pt,
  text(fill: white, size: 8pt, weight: "bold", label),
)

// KPI summary box used on the cover page: title, value, optional sub text
#let pcgr-kpi-box(title, value, sub: none, color: pcgr-neutral) = rect(
  fill: color,
  radius: 3pt,
  inset: (x: 10pt, y: 9pt),
  width: 100%,
  height: 2.3cm,
)[
  #set text(fill: white)
  // value (count, estimate or gene symbols) scaled down to fit the box width
  #layout(size => {
    let fs = 17pt
    let w = measure(text(size: fs, weight: "bold", value)).width
    if w > size.width { fs = calc.max(10pt, fs * (size.width / w)) }
    let items = (
      text(size: 8pt)[#title],
      v(12pt),
      text(size: fs, weight: "bold", value),
    )
    if sub != none and sub != "" {
      items.push(v(5pt))
      items.push(text(size: 8pt, sub))
    }
    stack(dir: ttb, ..items)
  })
]

// Metadata label+value pair (used inside the cover page info grid)
#let pcgr-meta(label, value) = (
  text(weight: "bold", size: 8.5pt, fill: pcgr-muted)[#upper(label)],
  text(size: 10pt, weight: "bold", fill: pcgr-dark)[#value],
)

#let pcgr-cover-heading(body) = text(
  size: 12pt, weight: "bold", fill: pcgr-main-dark)[#body]

// Placeholder for report content not yet implemented (skeleton)
#let pcgr-placeholder(body) = block(
  width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt,
  fill: pcgr-stripe, stroke: (paint: pcgr-rule, thickness: 0.8pt, dash: "dashed"),
  above: 0.8em, below: 1em,
)[
  #text(weight: "bold", fill: pcgr-main-dark)[To be implemented] #h(0.4em)
  #text(fill: pcgr-muted)[#body]
]

// ---------------------------------------------------------------------------
// Cover + executive summary page
// Called at the very top of pcgr_report_pdf.qmd via a raw Typst block emitted
// from R so that all values are resolved at render time.
// ---------------------------------------------------------------------------
#let pcgr-cover(
  sample_id:        "—",
  genome_assembly:  "—",
  assay:            "—",
  tumor_site:       "—",
  purity_ploidy:    "—",
  sex_age:          "—",
  pcgr_version:     "—",
) = {

  pcgr-sample.update(sample_id)
  pcgr-ver.update(pcgr_version)

  // --- Logo and report date ---
  grid(
    columns: (auto, 1fr),
    align: (left + horizon, right + horizon),
    pcgr-logo(2.4cm),
    align(right)[
      #set text(size: 8.5pt, fill: pcgr-muted)
      Report date
      #linebreak()
      #text(size: 10pt, weight: "bold", fill: pcgr-dark)[#pcgr-today]
    ],
  )

  v(0.6cm)

  // --- Full-bleed band: sample name as the headline ---
  pad(x: -pcgr-mx, rect(
    fill: pcgr-main, width: 100%, height: 5.2cm, stroke: none,
    inset: (x: pcgr-mx, bottom: 0.9cm, top: 0.9cm),
  )[
    #set text(fill: white)
    #align(bottom + left)[
      #text(size: 9pt, weight: "bold", tracking: 0.12em)[#upper[Cancer genome interpretation report]]
      #v(3pt)
      // sample ID on a single line: 26pt, scaled down to fit the band for
      // long IDs (min. 14pt)
      #layout(size => {
        let fs = 26pt
        let w = measure(text(size: fs, weight: "bold", sample_id)).width
        if w > size.width { fs = calc.max(14pt, fs * (size.width / w)) }
        text(size: fs, weight: "bold", sample_id)
      })
      #v(3pt)
      #text(size: 9pt)[Personal Cancer Genome Reporter · pcgr v#pcgr_version]
    ]
  ])

  v(0.8cm)

  // --- Sample metadata (two label/value columns) ---
  grid(
    columns: (auto, 1fr, auto, 1fr),
    column-gutter: 1.4em,
    row-gutter: 0.6em,
    ..pcgr-meta("Tumor site", tumor_site),
    ..pcgr-meta("Assay", assay),
    ..pcgr-meta("Sex / age", sex_age),
    ..pcgr-meta("Tumor purity / ploidy", purity_ploidy),
  )

  v(0.5cm)
  line(length: 100%, stroke: 0.5pt + pcgr-rule)
  v(0.4cm)
}

// Key findings on the cover: a framed panel with a colored title band, a
// tinted lead section (most important findings) and a white list of
// supporting findings. Content from pcgrr::report_highlights(); the R code
// emits #pcgr-highlights[ #pcgr-hl-lead[..markdown..] #pcgr-hl-support[..markdown..] ]
#let pcgr-hl-marker = (
  text(fill: pcgr-main, size: 8pt, weight: "bold")[#sym.bullet],
  text(fill: pcgr-main, weight: "bold")[–],
)
#let pcgr-highlights(body) = block(
  width: 100%, breakable: false, clip: true, radius: 5pt,
  stroke: 1pt + pcgr-main, fill: white,
)[
  #set block(spacing: 0pt)
  #block(width: 100%, fill: pcgr-main, inset: (x: 14pt, y: 8pt))[
    #set text(fill: white)
    #grid(
      columns: (auto, 1fr), column-gutter: 8pt, align: horizon,
      polygon.regular(vertices: 6, size: 11pt, fill: white),
      text(size: 10.5pt, weight: "bold", tracking: 0.14em)[#upper[Key findings]],
    )
  ]
  #body
]
#let pcgr-hl-lead(body) = block(
  width: 100%, inset: (x: 14pt, y: 12pt), fill: pcgr-main.lighten(90%),
)[
  #set text(size: 10.5pt)
  #set par(justify: false, leading: 0.62em)
  #set block(spacing: 0.8em)
  #set list(spacing: 0.65em, marker: pcgr-hl-marker)
  #show strong: set text(fill: pcgr-main-dark)
  #body
]
#let pcgr-hl-support(body) = block(
  width: 100%, inset: (x: 14pt, y: 11pt),
  stroke: (top: 0.75pt + pcgr-rule),
)[
  #set text(size: 9.5pt)
  #set par(justify: false, leading: 0.6em)
  #set block(spacing: 0.8em)
  #set list(spacing: 0.7em, marker: pcgr-hl-marker)
  #show strong: set text(fill: pcgr-main-dark)
  #body
]

// Positive composite genomic biomarkers (e.g. "TMB-high (12.3 mut/Mb)",
// "MSI-high") as pills, below the actionable alterations (lead section);
// only emitted by R when at least one biomarker is positive; each pill
// links to its report section when a target (heading label) is given
#let pcgr-hl-pills(labels, targets: ()) = block(above: 0.9em)[
  #for (i, l) in labels.enumerate() {
    // dark (not the oncogenicity/pathogenicity red) - positive composite
    // biomarkers may indicate treatment benefit
    let pill = box(fill: pcgr-dark, radius: 10pt, inset: (x: 11pt, y: 5pt))[
      #text(fill: white, size: 11pt, weight: "bold")[#l]
    ]
    if i < targets.len() { link(label(targets.at(i)), pill) } else { pill }
    h(7pt)
  }
]

// Actionable variants removed by variant filtering (depth/allelic support,
// tumor-only germline filtering), between the lead and the supporting
// findings; only emitted by R when such variants exist
#let pcgr-hl-warning(body) = block(
  width: 100%, inset: (x: 14pt, y: 10pt), fill: rgb("#fdf5e6"),
  stroke: (top: 0.75pt + pcgr-rule),
)[
  #set text(size: 9.5pt)
  #set par(justify: false, leading: 0.6em)
  #set block(spacing: 0.7em)
  #set list(spacing: 0.6em, marker: text(fill: rgb("#c87f0a"), weight: "bold")[!])
  #show strong: set text(fill: rgb("#8a5700"))
  #body
]

// Contents (level-1 headings, page numbers resolved by Typst) and end of
// the cover page
#let pcgr-cover-end() = {
  v(0.6cm)
  pcgr-cover-heading[Contents]
  v(0.2cm)
  {
    set text(size: 10pt)
    outline(title: none, depth: 1, indent: 0pt)
  }
  pagebreak()
}
