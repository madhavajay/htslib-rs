use htslib_rs::hfile_compat::{add_extension, decode_data_url};

#[test]
fn ports_haddextension_cases() {
    assert_eq!(
        add_extension("foo/bar.bam", false, ".bai"),
        "foo/bar.bam.bai"
    );
    assert_eq!(add_extension("foo/bar.bam", true, ".bai"), "foo/bar.bai");
    assert_eq!(
        add_extension("foo.bar/baz", true, ".bai"),
        "foo.bar/baz.bai"
    );
    assert_eq!(
        add_extension("foo#bar.bam", false, ".bai"),
        "foo#bar.bam.bai"
    );
    assert_eq!(add_extension(".bam", true, ".bai"), ".bai");
    assert_eq!(add_extension("foo", true, ".csi"), "foo.csi");
    assert_eq!(
        add_extension("http://host/bar.cram?a&b&c", false, ".crai"),
        "http://host/bar.cram.crai?a&b&c"
    );
    assert_eq!(
        add_extension("http://host/bar.cram#frag", true, ".crai"),
        "http://host/bar.crai#frag"
    );
}

#[test]
fn ports_data_url_hfile_cases() {
    assert_eq!(
        decode_data_url("data:,hello, world!%0A").unwrap(),
        b"hello, world!\n"
    );
    assert_eq!(decode_data_url("data:,").unwrap(), b"");
    assert_eq!(
        decode_data_url(concat!(
            "data:;base64,",
            "TWFuIGlzIGRpc3Rpbmd1aXNoZWQsIG5vdCBvbmx5IGJ5IGhpcyByZWFzb24sIGJ1dCBieSB0aGlz",
            "IHNpbmd1bGFyIHBhc3Npb24gZnJvbSBvdGhlciBhbmltYWxzLCB3aGljaCBpcyBhIGx1c3Qgb2Yg",
            "dGhlIG1pbmQsIHRoYXQgYnkgYSBwZXJzZXZlcmFuY2Ugb2YgZGVsaWdodCBpbiB0aGUgY29udGlu",
            "dWVkIGFuZCBpbmRlZmF0aWdhYmxlIGdlbmVyYXRpb24gb2Yga25vd2xlZGdlLCBleGNlZWRzIHRo",
            "ZSBzaG9ydCB2ZWhlbWVuY2Ugb2YgYW55IGNhcm5hbCBwbGVhc3VyZS4="
        ))
        .unwrap(),
        concat!(
            "Man is distinguished, not only by his reason, but by this singular passion from other ",
            "animals, which is a lust of the mind, that by a perseverance of delight in the continued ",
            "and indefatigable generation of knowledge, exceeds the short vehemence of any carnal pleasure."
        )
        .as_bytes()
    );
}
