use std::{
    io::{BufRead, Cursor, Read, Seek, SeekFrom, Write},
    num::NonZero,
};

use htslib_rs::{
    bgzf,
    bgzf_compat::{
        CompressionKind, build_gzi, detect_compression_kind, has_eof_marker, query_gzi, read_all,
        read_all_with_worker_count, read_auto, read_gzi, write_all, write_all_with_kind,
        write_all_with_worker_count, write_gzi,
    },
};

const PLAIN: &[u8] = include_bytes!("../../../repos/htslib/test/bgziptest.txt");
const BGZF: &[u8] = include_bytes!("../../../repos/htslib/test/bgziptest.txt.gz");
const GZI: &[u8] = include_bytes!("../../../repos/htslib/test/bgziptest.txt.gz.gzi");

fn generated_text() -> Vec<u8> {
    let mut text = Vec::new();

    for i in 0..50_000 {
        writeln!(&mut text, "{i:07}").expect("write to Vec");
    }

    text
}

#[test]
fn reads_upstream_bgzf_fixture() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(read_all(Cursor::new(BGZF))?, PLAIN);

    Ok(())
}

#[test]
fn writes_bgzf_with_eof_marker_and_reads_it_back() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let encoded = write_all(Vec::new(), &text)?;

    assert!(has_eof_marker(&encoded));
    assert_eq!(read_all(Cursor::new(encoded))?, text);

    Ok(())
}

#[test]
fn ports_bgzf_threaded_read_write_round_trip() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let worker_count = NonZero::new(2).unwrap();
    let encoded = write_all_with_worker_count(Vec::new(), &text, worker_count)?;

    assert!(has_eof_marker(&encoded));
    assert_eq!(
        read_all_with_worker_count(Cursor::new(encoded), worker_count)?,
        text
    );

    Ok(())
}

#[test]
fn ports_bgzf_plain_gzip_and_bgzf_write_modes() -> Result<(), Box<dyn std::error::Error>> {
    for kind in [
        CompressionKind::Uncompressed,
        CompressionKind::Gzip,
        CompressionKind::Bgzf,
    ] {
        let encoded = write_all_with_kind(PLAIN, kind)?;

        assert_eq!(detect_compression_kind(&encoded), kind);
        assert_eq!(read_auto(&encoded)?, PLAIN);
        assert_eq!(has_eof_marker(&encoded), kind == CompressionKind::Bgzf);
    }

    Ok(())
}

#[test]
fn reads_past_embedded_eof_block() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let mid = text.len() / 2;
    let mut encoded = write_all(Vec::new(), &text[..mid])?;
    encoded.extend(write_all(Vec::new(), &text[mid..])?);

    assert_eq!(encoded.windows(28).filter(|w| has_eof_marker(w)).count(), 2);
    assert_eq!(read_all(Cursor::new(encoded))?, text);

    Ok(())
}

#[test]
fn loads_dumps_and_queries_upstream_gzi() -> Result<(), Box<dyn std::error::Error>> {
    let index = read_gzi(Cursor::new(GZI))?;

    assert_eq!(write_gzi(Vec::new(), &index)?, GZI);

    for offset in [0, PLAIN.len() as u64 - 1] {
        let position = query_gzi(&index, offset)?;
        assert!(position.compressed() < BGZF.len() as u64);
    }

    Ok(())
}

#[test]
fn ports_test_rebgzip_multiblock_gzi_boundaries() -> Result<(), Box<dyn std::error::Error>> {
    let index = read_gzi(Cursor::new(GZI))?;

    assert_eq!(
        index.as_ref(),
        &[(29, 1), (59, 3), (90, 6), (122, 10), (153, 15)]
    );

    let mut reader = bgzf::io::IndexedReader::new(Cursor::new(BGZF), index.clone());

    for &(_, uncompressed_offset) in index.as_ref() {
        reader.seek(SeekFrom::Start(uncompressed_offset))?;

        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;

        assert_eq!(
            actual,
            PLAIN[usize::try_from(uncompressed_offset)?..],
            "uncompressed offset: {uncompressed_offset}"
        );
    }

    Ok(())
}

#[test]
fn seeks_by_uncompressed_offsets_with_gzi() -> Result<(), Box<dyn std::error::Error>> {
    let index = read_gzi(Cursor::new(GZI))?;
    let mut reader = bgzf::io::IndexedReader::new(Cursor::new(BGZF), index);

    for offset in [0, 5, 10] {
        reader.seek(SeekFrom::Start(offset as u64))?;

        let mut buf = [0; 5];
        reader.read_exact(&mut buf)?;

        assert_eq!(&buf, &PLAIN[offset..offset + buf.len()]);
    }

    Ok(())
}

#[test]
fn seeks_by_virtual_positions_recorded_while_writing() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let chunk_len = text.len() / 10;
    let mut writer = bgzf::io::Writer::new(Vec::new());
    let mut positions = Vec::new();

    for chunk in text.chunks(chunk_len).take(10) {
        positions.push(writer.virtual_position());
        writer.write_all(chunk)?;
    }

    let encoded = writer.finish()?;
    let mut reader = bgzf::io::Reader::new(Cursor::new(encoded));

    for (i, position) in positions.into_iter().enumerate() {
        reader.seek(position)?;

        let mut buf = vec![0; 16];
        reader.read_exact(&mut buf)?;

        let offset = i * chunk_len;
        assert_eq!(buf, text[offset..offset + 16]);
    }

    Ok(())
}

#[test]
fn reads_lines_like_bgzf_getline() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let encoded = write_all(Vec::new(), &text)?;
    let mut reader = bgzf::io::Reader::new(Cursor::new(encoded));

    for line in text.split(|b| *b == b'\n').filter(|line| !line.is_empty()) {
        let mut buf = Vec::new();
        let n = reader.read_until(b'\n', &mut buf)?;

        assert_eq!(n, line.len() + 1);
        assert_eq!(buf.strip_suffix(b"\n").unwrap_or(&buf), line);
    }

    assert!(reader.fill_buf()?.is_empty());

    Ok(())
}

#[test]
fn rejects_truncated_bgzf_stream() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let mut encoded = write_all(Vec::new(), &text)?;
    encoded.truncate(encoded.len() / 2);

    assert!(read_all(Cursor::new(encoded)).is_err());

    Ok(())
}

#[test]
fn builds_gzi_for_written_bgzf() -> Result<(), Box<dyn std::error::Error>> {
    let text = generated_text();
    let encoded = write_all(Vec::new(), &text)?;
    let index = build_gzi(&mut Cursor::new(&encoded))?;
    let mut reader = bgzf::io::IndexedReader::new(Cursor::new(encoded), index.clone());

    for offset in [0, 100, 50, 70_000, text.len() - 16] {
        reader.seek(SeekFrom::Start(offset as u64))?;

        let mut buf = [0; 16];
        reader.read_exact(&mut buf)?;

        assert_eq!(&buf, &text[offset..offset + buf.len()]);
    }

    assert!(query_gzi(&index, 70_000)?.compressed() > 0);

    Ok(())
}
