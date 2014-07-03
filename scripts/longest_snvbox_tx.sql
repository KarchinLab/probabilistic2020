/* Retrieves the longest transcript for each gene in SNVBox.

The longest RefSeq transcript is selected first. However, if there are no
RefSeq transcripts for a gene, then the longest Ensembl transcript is selected.

To output results to a file use mysql from the command line.

$ mysql [options] < longest_snvbox_tx.sql > output.txt

Please change the database name your_snvbox_db to your actual SNVBox db name.
*/
USE SNVBox_dev;  -- use your own db name for SNVBox
SELECT * 
FROM (
    SELECT c.GeneSymbol as `Gene Symbol`, c.RefseqT as `Transcript`, c.aaLen as `AA Length`, b.max_aaLen as `Max AA Length`
    FROM (
        SELECT c_gs1.UID, c_gs1.GeneSymbol, c_t1.RefseqT, c_t1.EnsT, c_t1.aaLen
        FROM GeneSymbols c_gs1
        INNER JOIN (
            SELECT *
            FROM Transcript
            WHERE Transcript.RefseqT IS NOT NULL
        ) c_t1
        ON c_t1.UID=c_gs1.UID
    ) c
    INNER JOIN (
        SELECT a.GeneSymbol, MAX(a.aaLen) as max_aaLen
        FROM (
            SELECT a_gs2.UID, a_gs2.GeneSymbol, a_t2.RefseqT, a_t2.aaLen
            FROM GeneSymbols a_gs2 
            INNER JOIN ( 
                SELECT *
                FROM Transcript
                WHERE Transcript.RefseqT IS NOT NULL 
            ) a_t2
            ON a_gs2.UID=a_t2.UID
        ) a 
        GROUP BY a.GeneSymbol
    ) b
    ON b.GeneSymbol=c.GeneSymbol AND b.max_aaLen=c.aaLen
) lr_only
UNION (
    SELECT le.*
    FROM (
        SELECT z.GeneSymbol as `Gene Symbol`, z.EnsT as `Transcript`, z.aaLen as `AA Length`, y.max_aaLen as `Max AA Length`
        FROM (
            SELECT gs1.UID, gs1.GeneSymbol, t1.RefseqT, t1.EnsT, t1.aaLen
            FROM GeneSymbols gs1
            INNER JOIN (
                SELECT *
                FROM Transcript
                WHERE Transcript.EnsT IS NOT NULL
            ) t1
            ON t1.UID=gs1.UID
        ) z
        INNER JOIN (
            SELECT x.GeneSymbol, MAX(x.aaLen) as max_aaLen
            FROM (
                SELECT gs2.UID, gs2.GeneSymbol, t2.EnsT, t2.aaLen
                FROM GeneSymbols gs2 
                INNER JOIN Transcript t2
                ON gs2.UID=t2.UID
            ) x
            GROUP BY x.GeneSymbol
        ) y
        ON y.GeneSymbol=z.GeneSymbol AND y.max_aaLen=z.aaLen
    ) le
    LEFT JOIN (
        SELECT zz.GeneSymbol as `Gene Symbol`, zz.RefseqT as `Transcript`, zz.aaLen as `AA Length`, yy.max_aaLen as `Max AA Length`
        FROM (
            SELECT zz_gs1.UID, zz_gs1.GeneSymbol, zz_t1.RefseqT, zz_t1.EnsT, zz_t1.aaLen
            FROM GeneSymbols zz_gs1
            INNER JOIN (
                SELECT *
                FROM Transcript
                WHERE Transcript.RefseqT IS NOT NULL
            ) zz_t1
            ON zz_t1.UID=zz_gs1.UID
        ) zz
        INNER JOIN (
            SELECT xx.GeneSymbol, MAX(xx.aaLen) as max_aaLen
            FROM (
                SELECT xx_gs2.UID, xx_gs2.GeneSymbol, xx_t2.RefseqT, xx_t2.aaLen
                FROM GeneSymbols xx_gs2 
                INNER JOIN (
                    SELECT *
                    FROM Transcript 
                    WHERE Transcript.RefseqT IS NOT NULL
                ) xx_t2
                ON xx_gs2.UID=xx_t2.UID
            ) xx
            GROUP BY xx.GeneSymbol
        ) yy
        ON yy.GeneSymbol=zz.GeneSymbol AND yy.max_aaLen=zz.aaLen
    ) lr
    ON lr.`Gene Symbol`=le.`Gene Symbol`
    WHERE lr.Transcript IS NULL
)
ORDER BY `Gene Symbol`;
