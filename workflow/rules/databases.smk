# The download rules exist only when db_dir is "resources/db", so a run never writes to a shared
# database folder. A shared folder is filled once with the download_databases target.
if not DB_SHARED:

  localrules: download_databases

  rule download_sprot:
    """
    Downloads the Swiss-Prot release named by sprot_release, checks its md5 and prepares it for BLAST.
    UniProt serves only the newest release as a plain fasta file. An older release comes from the
    release archive, which holds the same fasta file.
    """
    output:
      fasta="resources/db/uniprot_sprot.fasta",
      version="resources/db/uniprot_sprot.version"
    log:
      "logs/download_sprot.log"
    conda:
      "../envs/blast.yaml"
    params:
      release=SPROT_RELEASE,
      current_url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete",
      archive_url="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-" + SPROT_RELEASE + "/knowledgebase"
    threads:
      1
    resources:
      mem_mb=2048,
      runtime="4h"
    shell:
      """
      db=$(dirname {output.fasta})

      # the release that UniProt serves as current
      current_release() {{
        wget -q -O - "{params.current_url}/reldate.txt" | awk '$1 == "UniProtKB/Swiss-Prot" && $2 == "Release" {{ print $3 }}'
      }}

      # the md5 that the RELEASE.metalink file in folder $1 gives for file $2
      metalink_md5() {{
        wget -q -O - "$1/RELEASE.metalink" | awk -v name="$2" '/<file name=/ {{ f = index($0, "\\"" name "\\"") > 0 }} f && !found && /type="md5"/ {{ split($0, a, /[<>]/); print a[3]; found = 1 }}'
      }}

      before=$(current_release)
      echo "Swiss-Prot release asked: {params.release}, current release: $before" > {log}

      if [ "$before" = "{params.release}" ]; then
        download="$db/uniprot_sprot.fasta.gz"
        wget -O "$download" "{params.current_url}/uniprot_sprot.fasta.gz" &>> {log}
        expected=$(metalink_md5 "{params.current_url}" uniprot_sprot.fasta.gz)
      else
        archive="uniprot_sprot-only{params.release}.tar.gz"
        download="$db/$archive"
        if ! wget -O "$download" "{params.archive_url}/$archive" &>> {log}; then
          echo "Swiss-Prot {params.release} is not the current release ($before) and has no archive at {params.archive_url}" >> {log}
          exit 1
        fi
        expected=$(metalink_md5 "{params.archive_url}" "$archive")
      fi

      actual=$(md5sum "$download" | cut -d " " -f 1)
      echo "md5 of $download: expected $expected, actual $actual" >> {log}
      if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
        echo "The md5 does not match RELEASE.metalink" >> {log}
        exit 1
      fi

      if [ "$before" = "{params.release}" ]; then
        # UniProt could switch to the next release while the file downloads
        after=$(current_release)
        if [ "$after" != "$before" ]; then
          echo "The current release changed from $before to $after during the download. Run the rule again." >> {log}
          exit 1
        fi
      else
        tar -xzf "$download" -C "$db" uniprot_sprot.fasta.gz &>> {log}
        rm "$download"
      fi

      gunzip -f "$db/uniprot_sprot.fasta.gz" &>> {log}
      makeblastdb -in {output.fasta} -dbtype prot &>> {log}
      printf "Swiss-Prot\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_pfam:
    """
    Downloads the Pfam release named by pfam_release, checks its md5 and prepares it for hmmscan.
    """
    output:
      hmm="resources/db/Pfam-A.hmm",
      version="resources/db/Pfam-A.version"
    log:
      "logs/download_pfam.log"
    conda:
      "../envs/pfam.yaml"
    params:
      release=PFAM_RELEASE,
      url="https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam" + PFAM_RELEASE
    threads:
      1
    resources:
      mem_mb=2048,
      runtime="4h"
    shell:
      """
      wget -O {output.hmm}.gz "{params.url}/Pfam-A.hmm.gz" &> {log}
      expected=$(wget -q -O - "{params.url}/md5_checksums" | awk '$2 == "Pfam-A.hmm.gz" {{ print $1 }}')
      actual=$(md5sum {output.hmm}.gz | cut -d " " -f 1)
      echo "md5 of Pfam-A.hmm.gz: expected $expected, actual $actual" >> {log}
      if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
        echo "The md5 does not match md5_checksums of Pfam {params.release}" >> {log}
        exit 1
      fi
      gunzip -f {output.hmm}.gz &>> {log}
      # -f replaces index files left by an earlier attempt
      hmmpress -f {output.hmm} &>> {log}
      printf "Pfam\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_rfam:
    """
    Downloads the Rfam release named by rfam_release and prepares it for cmscan.
    Rfam publishes no checksum, so the file is not checked.
    """
    output:
      cm="resources/db/Rfam.cm",
      version="resources/db/Rfam.version"
    log:
      "logs/download_rfam.log"
    conda:
      "../envs/rfam.yaml"
    params:
      release=RFAM_RELEASE,
      url="https://ftp.ebi.ac.uk/pub/databases/Rfam/" + RFAM_RELEASE
    threads:
      1
    resources:
      mem_mb=2048,
      runtime="4h"
    shell:
      """
      wget -O {output.cm}.gz "{params.url}/Rfam.cm.gz" &> {log}
      gunzip -f {output.cm}.gz &>> {log}
      # -F replaces index files left by an earlier attempt
      cmpress -F {output.cm} &>> {log}
      printf "Rfam\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_busco_lineage:
    """
    Downloads the BUSCO lineage named by lineage into resources/db/busco_downloads, for a shared database
    folder. A run with db_dir "resources/db" does not read this copy: BUSCO downloads its own lineage
    into resources/busco_downloads.
    """
    output:
      directory("resources/db/busco_downloads/lineages/" + LINEAGE)
    log:
      "logs/download_busco_lineage.log"
    conda:
      "../envs/busco.yaml"
    params:
      lineage=LINEAGE
    threads:
      1
    resources:
      mem_mb=2048,
      runtime="4h"
    shell:
      """
      # --download has to come last on the busco command line
      busco --download_path resources/db/busco_downloads --download {params.lineage} &> {log}
      """


  rule download_databases:
    """
    Downloads the annotation databases and the BUSCO lineage into resources/db/ and lists their releases
    in resources/db/versions.txt. Run it once in an empty folder to fill a shared database folder.
    """
    input:
      versions=["resources/db/uniprot_sprot.version", "resources/db/Pfam-A.version", "resources/db/Rfam.version"],
      lineage=["resources/db/busco_downloads/lineages/" + LINEAGE] if LINEAGE else []
    output:
      "resources/db/versions.txt"
    threads:
      1
    resources:
      mem_mb=2048,
      runtime="4h"
    run:
      with open(output[0], "w") as out:
        for path in input.versions:
          with open(path) as handle:
            out.write(handle.read())
        for path in input.lineage:
          fields = {}
          with open(os.path.join(path, "dataset.cfg")) as handle:
            for line in handle:
              if "=" in line:
                key, value = line.strip().split("=", 1)
                fields[key] = value
          out.write("BUSCO\t%s\t%s\n" % (fields.get("name", os.path.basename(path)), fields.get("creation_date", "unknown")))
