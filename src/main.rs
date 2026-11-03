//! Group BAM records into read groups by the value of a tag.
use std::num::NonZero;
use std::path::PathBuf;
use std::process;

use anyhow::Error;
use bamby::{Args, run};
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::{CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use log::error;
use mimalloc::MiMalloc;
use noodles::sam::alignment::record::data::field::Tag;
use noodles_bgzf::io::writer::CompressionLevel;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// The tool's one-line description at the top of the help: magenta, the same
/// color code as the [`FOOTER`], so the help opens and closes in one color.
pub(crate) const TITLE: Style = AnsiColor::Magenta.on_default();
/// Indented code examples and backtick-wrapped terms in help text.
pub(crate) const CODE: Style = AnsiColor::Yellow.on_default();
/// The license footer at the bottom of the help, shared with the [`TITLE`].
pub(crate) const FOOTER: Style = AnsiColor::Magenta.on_default();
/// The right-hand description column of an indented two-column table.
pub(crate) const FADED: Style = AnsiColor::BrightBlack.on_default();

/// Cargo's color style.
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Group BAM records into read groups by the value of a tag.
///
/// bamby gives every distinct value of a tag, such as the cell barcode in
/// `CB`, a read group of its own, so tools that work per read group or per
/// sample see one sample per value: one BAM holds every cell, and
/// `bcftools mpileup` reports each cell as its own sample. The input is read
/// twice: once to learn every value and build the header, then again to write
/// each record with its new `RG` tag.
///
/// READ GROUPS
///
///   ID    the record's read group and the value, joined by a period
///         (`A.AACGTT`), or the value alone when the record has none
///   SM    the value
///   rest  every other field (`LB`, `PL`, `PU`, ...) copied from the
///         record's read group
///
/// Records without the tag keep their read group, or are dropped with
/// `--drop-untagged`. A read group no written record uses any more is removed
/// from the header, and new read groups are listed in sorted order.
///
/// EXAMPLES
///
///   bamby -i in.bam -o out.bam
///   bamby -i in.bam -o - -t CB -d -l 0 | samtools sort -o sorted.bam
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    color = clap::ColorChoice::Always,
    rename_all = "kebab-case",
    verbatim_doc_comment,
    arg_required_else_help = true
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    /// Input BAM, in any order.
    ///
    /// Read twice, so it must be a file.
    #[arg(short = 'i', long, value_name = "BAM", value_parser = file, verbatim_doc_comment)]
    input: PathBuf,

    /// Output BAM, in the input's order; `-` writes to standard output.
    #[arg(short = 'o', long, value_name = "BAM", verbatim_doc_comment)]
    output: PathBuf,

    /// The tag whose values name the read groups.
    #[arg(
        short = 't',
        long,
        value_name = "TAG",
        default_value = "CB",
        value_parser = tag,
        verbatim_doc_comment
    )]
    tag: Tag,

    /// Drop records without the tag instead of keeping their read group.
    #[arg(short = 'd', long, verbatim_doc_comment)]
    drop_untagged: bool,

    /// BGZF worker threads, for reading and again for writing.
    #[arg(
        short = '@',
        long,
        value_name = "N",
        default_value = "1",
        verbatim_doc_comment
    )]
    threads: NonZero<usize>,

    /// BGZF compression level of the output, from 0 (none) to 12 (best).
    ///
    /// Use 0 when piping into another tool such as `samtools sort`.
    #[arg(
        short = 'l',
        long,
        value_name = "LEVEL",
        default_value = "6",
        value_parser = compression_level,
        verbatim_doc_comment
    )]
    compression_level: CompressionLevel,
}

fn file(text: &str) -> Result<PathBuf, String> {
    match text {
        "-" => Err("the input is read twice, so it must be a file, not standard input".into()),
        _ => Ok(PathBuf::from(text)),
    }
}

fn tag(text: &str) -> Result<Tag, String> {
    let bytes = text.as_bytes();
    let valid =
        bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric();
    if valid {
        Ok(Tag::new(bytes[0], bytes[1]))
    } else {
        Err(format!(
            "a SAM tag is a letter then a letter or digit, like CB, not {text:?}"
        ))
    }
}

fn compression_level(text: &str) -> Result<CompressionLevel, String> {
    text.parse::<u8>()
        .ok()
        .and_then(CompressionLevel::new)
        .ok_or_else(|| format!("a compression level is a whole number from 0 to 12, not {text:?}"))
}

/// The ANSI escape that starts `style`, or an empty string when `color` is off
/// (honoring `NO_COLOR`). Its matching reset comes from [`esc_reset`].
fn esc(style: Style, color: bool) -> String {
    if color {
        style.render().to_string()
    } else {
        String::new()
    }
}

/// The reset escape for `style` (a full SGR reset), or empty when `color` is off.
fn esc_reset(style: Style, color: bool) -> String {
    if color {
        style.render_reset().to_string()
    } else {
        String::new()
    }
}

/// Paint the first line of `text` in the [`TITLE`] color (the tool's one-line
/// description), leaving the rest of the text untouched.
fn title_first_line(text: &str, color: bool) -> String {
    let title = esc(TITLE, color);
    let reset = esc_reset(TITLE, color);
    match text.split_once('\n') {
        Some((first, rest)) => format!("{title}{first}{reset}\n{rest}"),
        None => format!("{title}{text}{reset}"),
    }
}

/// Drop the backticks around terms (`` `like this` ``), painting each term with
/// the `code` escape. `after` is the escape that restores the surrounding text's
/// color once a term ends.
fn paint_backtick_terms(text: &str, code: &str, after: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        out.push_str(&rest[..open]);
        let tail = &rest[open + 1..];
        if let Some(close) = tail.find('`') {
            out.push_str(code);
            out.push_str(&tail[..close]);
            out.push_str(after);
            rest = &tail[close + 1..];
        } else {
            out.push('`');
            rest = tail;
            break;
        }
    }
    out.push_str(rest);
    out
}

/// Split an indented two-column row into (left, gap, right) at the first run of
/// two or more spaces that follows the left-hand term.
fn split_two_column(content: &str) -> Option<(&str, &str, &str)> {
    let bytes = content.as_bytes();
    let mut seen_term = false;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b' ' {
            let start = i;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            if seen_term && i - start >= 2 {
                return Some((&content[..start], &content[start..i], &content[i..]));
            }
        } else {
            seen_term = true;
            i += 1;
        }
    }
    None
}

/// Style one help string, dropping backticks throughout:
///
/// - indented two-column rows: the left code term in [`CODE`], the right
///   description column in [`FADED`];
/// - indented description continuations (deeper indent, no term): all [`FADED`];
/// - indented single-column examples (a bare command): all [`CODE`];
/// - prose: default color, with backtick terms painted in [`CODE`].
fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            if !content.starts_with("  ") || content.trim().is_empty() {
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            let indent = content.len() - content.trim_start().len();
            if indent > 2 {
                let painted = paint_backtick_terms(content, &code, &faded);
                format!("{faded}{painted}{reset}{newline}")
            } else if content.trim_start().starts_with(env!("CARGO_PKG_NAME")) {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            } else if let Some((left, gap, right)) = split_two_column(content) {
                let left = paint_backtick_terms(left, &code, &code);
                let right = paint_backtick_terms(right, &code, &faded);
                format!("{code}{left}{reset}{gap}{faded}{right}{reset}{newline}")
            } else {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            }
        })
        .collect()
}

/// Style a command's help: paint the first line of the about text in the
/// [`TITLE`] color, apply [`style_help_text`] to the abouts and every option's
/// help, and add the license footer.
fn decorate_help(cmd: clap::Command, color: bool) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        esc(FOOTER, color),
        esc_reset(FOOTER, color),
    );

    let choice = if color {
        clap::ColorChoice::Always
    } else {
        clap::ColorChoice::Never
    };
    let mut cmd = cmd.color(choice);
    if let Some(about) = about {
        cmd = cmd.about(title_first_line(&style_help_text(&about, color), color));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(title_first_line(
            &style_help_text(&long_about, color),
            color,
        ));
    }
    cmd = cmd
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help, color));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help, color));
        }
        arg
    })
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let color = std::env::var_os("NO_COLOR").is_none();

    let env = Env::default().default_filter_or("info");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Help text is hand-wrapped, and clap would count the injected ANSI escapes as width.
    let mut cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches_mut();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    let args = Args {
        input: cli.input,
        output: cli.output,
        tag: cli.tag,
        drop_untagged: cli.drop_untagged,
        threads: cli.threads,
        compression_level: cli.compression_level,
        command_line: std::env::args().collect::<Vec<_>>().join(" "),
    };

    match run(&args) {
        Ok(_) => process::exit(0),
        Err(e) => {
            error!("{e:#}");
            process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use rstest::rstest;

    #[rstest]
    #[case("CB")]
    #[case("RX")]
    #[case("X0")]
    fn test_tag_accepts_valid_sam_tags(#[case] text: &str) {
        assert!(tag(text).is_ok());
    }

    #[rstest]
    #[case("C")]
    #[case("CBX")]
    #[case("0B")]
    #[case("C-")]
    fn test_tag_rejects_invalid_sam_tags(#[case] text: &str) {
        assert!(tag(text).is_err());
    }

    #[rstest]
    #[case("0", Some(0))]
    #[case("12", Some(12))]
    #[case("13", None)]
    #[case("fast", None)]
    fn test_compression_level_parses_zero_through_twelve(
        #[case] text: &str,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(
            compression_level(text).ok().map(|level| level.get()),
            expected
        );
    }

    #[test]
    fn test_file_rejects_standard_input() {
        assert!(file("-").is_err());
        assert!(file("in.bam").is_ok());
    }

    #[test]
    fn test_two_column_rows_split_at_the_first_wide_gap() {
        assert_eq!(
            split_two_column("  ID    the value"),
            Some(("  ID", "    ", "the value"))
        );
        assert_eq!(split_two_column("  bamby -i in.bam"), None);
    }

    #[test]
    fn test_style_help_text_is_plain_without_color() {
        let text = "Use `CB`.\n  ID    the value\n";
        assert_eq!(style_help_text(text, false), "Use CB.\n  ID    the value\n");
    }
}
