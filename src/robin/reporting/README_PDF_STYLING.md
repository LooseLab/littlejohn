# ROBIN PDF Report Styling System

## Overview

The ROBIN PDF report styling system has been completely updated to implement **Material Design 3 (M3)** and **Apple Human Interface Guidelines (HIG)** principles. This creates modern, professional, and accessible PDF reports that maintain consistency with the web interface while ensuring optimal readability and visual hierarchy.

## Design Principles

### Material Design 3 Integration
- **Dynamic Color System**: Uses semantic color tokens for consistent theming
- **Typography Scale**: Implements M3 typography hierarchy (Display, Headline, Title, Body, Label)
- **Spacing System**: 8dp grid-based spacing for consistent layouts
- **Elevation**: Subtle shadows and layering for depth perception
- **Component Styling**: Modern table designs, cards, and interactive elements

### Apple HIG Integration
- **Clarity**: High contrast and readable typography
- **Deference**: Content-focused design with subtle visual elements
- **Depth**: Layered information architecture
- **Spring Physics**: Natural, responsive animations and transitions

## Color System

### Primary Colors (ROBIN Green Theme)
```python
primary: "#4F9153"           # Main green
primary_container: "#E8F5E8"  # Light green background
on_primary: "#FFFFFF"        # White text on green
```

### Secondary Colors
```python
secondary: "#6750A4"         # M3 secondary purple
secondary_container: "#EADDFF" # Light purple background
on_secondary: "#FFFFFF"      # White text on purple
```

### Surface Colors
```python
surface: "#FFFFFF"           # Pure white background
surface_variant: "#F3F4F6"   # Light gray background
on_surface: "#1C1B1F"       # Dark text on white
on_surface_variant: "#49454F" # Medium text on white
```

### Semantic Colors
```python
success: "#2E7D32"           # Success green
warning: "#F57C00"           # Warning orange
error: "#C62828"             # Error red
info: "#1976D2"              # Info blue
```

## Typography System

### Display Styles
- **DisplayLarge**: 22pt, bold, centered - Main report titles
- **DisplayMedium**: 18pt, bold, centered - Section headers
- **DisplaySmall**: 16pt, bold, centered - Subsection headers

### Headline Styles
- **HeadlineLarge**: 16pt, bold - Major section headers
- **HeadlineMedium**: 14pt, medium - Section headers
- **HeadlineSmall**: 12pt, medium - Subsection headers

### Title Styles
- **TitleLarge**: 13pt, medium - Card titles, table headers
- **TitleMedium**: 12pt, medium - Subsection titles
- **TitleSmall**: 11pt, medium - Minor titles

### Body Styles
- **BodyLarge**: 11pt, regular - Main content text
- **BodyMedium**: 10pt, regular - Secondary content
- **BodySmall**: 9pt, regular - Supporting text

### Label Styles
- **LabelLarge**: 10pt, medium - Form labels, captions
- **LabelMedium**: 9pt, medium - Secondary labels
- **LabelSmall**: 8pt, medium - Fine print, metadata

## Component Styling

### Tables
Three table styles are available:

1. **MODERN_TABLE_STYLE**: Standard tables with alternating row colors
2. **COMPACT_TABLE_STYLE**: Dense tables for data-heavy content
3. **SUMMARY_TABLE_STYLE**: Highlighted tables for key information

### Cards and Containers
- **SummaryCard**: Information cards with borders and backgrounds
- **Alert**: Warning and information boxes
- **Success/Warning/Error**: Status-specific styling

### Code and Monospace
- **MonoText**: Inline code and technical data
- **CodeBlock**: Code blocks with language support

## Usage Examples

### Basic Text Styling
```python
from robin.reporting.styling import PDFStyleUtils

# Create section headers
header = PDFStyleUtils.create_section_header("Analysis Results", "HeadlineLarge")

# Create body text
body = PDFStyleUtils.create_body_text("This is the main content.")

# Create captions
caption = PDFStyleUtils.create_caption("Figure 1: Sample distribution")
```

### Table Creation
```python
# Create data table with headers
headers = ["Sample", "Quality", "Coverage"]
data = [["Sample1", "High", "95%"], ["Sample2", "Medium", "87%"]]
table = PDFStyleUtils.create_data_table(data, headers, "MODERN_TABLE_STYLE")
```

### Information Cards
```python
# Create info card
content = ["High quality reads", "Good coverage", "No contamination detected"]
card = PDFStyleUtils.create_info_card("Quality Summary", content, "SummaryCard")
```

### Status Indicators
```python
# Create status indicators
success_msg = PDFStyleUtils.create_status_indicator("success", "Analysis completed successfully")
warning_msg = PDFStyleUtils.create_status_indicator("warning", "Low coverage detected")
```

### Layout Utilities
```python
from robin.reporting.styling import PDFLayoutUtils

# Create two-column layout
left_col = [header1, content1]
right_col = [header2, content2]
layout = PDFLayoutUtils.create_two_column_layout(left_col, right_col)
```

## Font System

### Primary Fonts
- **Inter**: Modern, readable sans-serif for body text
- **JetBrains Mono**: Monospace font for code and technical data

### Fallback Fonts
- **Helvetica**: System fallback for body text
- **Courier**: System fallback for monospace

## Spacing System

### M3 Spacing Units
- **8pt**: Small spacing between related elements
- **16pt**: Standard spacing between sections
- **24pt**: Large spacing between major sections
- **32pt**: Extra large spacing for page breaks

### Implementation
```python
from reportlab.platypus import Spacer

# Standard spacing
spacer = Spacer(1, 16)  # 16pt spacing

# Small spacing
small_spacer = Spacer(1, 8)  # 8pt spacing
```

## Best Practices

### Typography
1. Use appropriate heading hierarchy (Display → Headline → Title → Body)
2. Maintain consistent spacing between related elements
3. Ensure sufficient contrast for readability

### Layout
1. Use the 8pt grid system for consistent spacing
2. Group related information in cards or containers
3. Use tables for structured data display

### Color Usage
1. Use semantic colors consistently (success, warning, error)
2. Maintain high contrast ratios for accessibility
3. Use primary colors sparingly for emphasis

### Content Organization
1. Start with clear section headers
2. Use cards for important information
3. Include appropriate spacing between sections
4. Use status indicators for key results

## Migration Guide

### From Old Styling
1. Replace old style names with new M3 equivalents:
   - `"Title"` → `"DisplayLarge"`
   - `"Heading1"` → `"HeadlineLarge"`
   - `"Heading2"` → `"HeadlineMedium"`
   - `"Normal"` → `"BodyLarge"`

2. Update table styling:
   - Use `MODERN_TABLE_STYLE` for standard tables
   - Use `COMPACT_TABLE_STYLE` for dense data
   - Use `SUMMARY_TABLE_STYLE` for key information

3. Implement consistent spacing:
   - Use `PDFStyleUtils.create_section_spacing()` between sections
   - Use `Spacer(1, 16)` for standard spacing
   - Use `Spacer(1, 8)` for small spacing

### New Features
1. **Utility Classes**: Use `PDFStyleUtils` for common styling tasks
2. **Layout Helpers**: Use `PDFLayoutUtils` for complex layouts
3. **Color Management**: Use `PDFColorUtils` for consistent color usage
4. **Status Indicators**: Use semantic styling for different message types

## Accessibility Features

- **High Contrast**: Ensures text readability
- **Semantic Colors**: Consistent meaning across the interface
- **Clear Typography**: Readable font sizes and weights
- **Logical Structure**: Clear heading hierarchy
- **Consistent Spacing**: Predictable layout patterns

## Performance Considerations

- Font loading is optimized with fallbacks
- Styles are cached for efficient reuse
- Minimal overhead for styling operations
- Efficient table rendering with pre-defined styles

## Future Enhancements

- Dark mode support
- Custom color scheme support
- Additional component styles
- Interactive PDF elements
- Enhanced accessibility features

---

For questions or contributions to the styling system, please refer to the main ROBIN documentation or contact the development team.
