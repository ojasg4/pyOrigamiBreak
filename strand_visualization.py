"""
strand_visualization.py

Interactive strand-level visualization tool for DNA origami designs.
This module creates interactive heatmaps where individual strands (sections of oligos on single helices)
can be clicked to display their thermodynamic properties.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
import xml.etree.ElementTree as ET
import re
import html
from matplotlib.colors import to_hex


class StrandVisualizer:
    """
    Creates strand-level interactive visualizations for DNA origami designs.
    """
    
    def __init__(self, autobreak_instance, override):
        """
        Initialize with an AutoBreak instance that has completed processing.
        
        Args:
            autobreak_instance: AutoBreak object with processed origami data
        """
        self.autobreak = autobreak_instance
        self.origami = autobreak_instance.origami
        if override is None:
            self.output_directory = autobreak_instance.output_directory
            self.output_name = autobreak_instance.output_name
            self.strand_data_file = os.path.join(self.output_directory, 'outputs', 
                                           f'{self.output_name}_strand_data.xlsx')
            self.strand_svg_path = os.path.join(self.output_directory, 'intermediates',
                                            f'{self.output_name}_strand_heatmap.svg')
            self.interactive_html_path = os.path.join(self.output_directory, 'outputs',
                                                    f'{self.output_name}_interactive_strands.html')
        else:
            self.output_directory = override[0]
            self.output_name = override[1]
            self.strand_data_file = os.path.join(self.output_directory, 'inputs', 'original_strand_data.xlsx')
            self.strand_svg_path = os.path.join(self.output_directory, 'inputs', 'original_strand_heatmap.svg')
            self.interactive_html_path = os.path.join(self.output_directory, 'inputs', 'original_interactive_strands.html')
    
    def extract_strand_data(self):
        """
        Extract thermodynamic data for individual strands from the origami structure.
        Returns a list of dictionaries containing strand-level data.
        """
        strand_data = []
        strand_id = 0
        
        for oligo in self.origami.oligos['staple']:
            current_strand = oligo.null_strand.next_strand
            
            while current_strand:
                # Calculate strand-level thermodynamics
                strand_info = self._calculate_strand_thermodynamics(current_strand, oligo, strand_id)
                if strand_info:
                    strand_data.append(strand_info)
                    strand_id += 1
                
                current_strand = current_strand.next_strand
        min, max = self._calculate_min_max_tf(strand_data)
        for strand in strand_data:
            strand['TfColor'] = self._get_tf_color(strand['Tf'], min, max)
        # strand_data.sort(key=lambda x: x['Tf'])
        # for i, strand in enumerate(strand_data):
        #     strand['StrandID'] = i
        # FIX
        return strand_data
    
    def _calculate_min_max_tf(self, strand_data):
        "Helper function to generate the max and min bounds for relative color assignment. +-1 to have open boundaries"
        min = 50.0
        max = 50.0
        for strand in strand_data:
            tf = strand['Tf']
            if tf < min:
                min = tf
            if tf > max:
                max = tf
        max += 1
        min -= 1
        return min, max

    def _calculate_strand_thermodynamics(self, strand, oligo, strand_id):
        """
        Calculate thermodynamic properties for a single strand using utilities directly.
        
        Args:
            strand: Strand object
            oligo: Parent oligo object
            strand_id: Unique identifier for this strand
            
        Returns:
            Dictionary with strand thermodynamic data
        """
        if not hasattr(strand, 'dna') or not strand.dna:
            print(f"DEBUG: Strand {strand_id} has no DNA sequence")
            return None
        
        # Get strand sequence and clean it
        sequence = strand.dna.replace(' ', '').replace('?', '')
        if not sequence:
            return None
        
        try:
            from autobreak import utilities
            
            # Get temperature from autobreak settings
            temperature_kelvin = self.autobreak.optim_temperature_kelvin
            
            # Calculate basic thermodynamics using utilities functions
            dG_total, dH_total, dS_total = utilities.sequence_to_dG_dH_dS(sequence, temperature_kelvin)
            tm_value = utilities.sequence_to_Tm(sequence)
            dG_conc, dS_conc = utilities.conc_to_dG(temperature_kelvin)
            
            # Calculate folding probability and Tf
            RT = utilities.R * temperature_kelvin
            edge_prob = np.exp(-dG_total/RT) / (1.0 + np.exp(-dG_total/RT))
            edge_logprob = np.log(edge_prob) if edge_prob > 0 else -50
            edge_Tf = dH_total/dS_total - 273.15 if dS_total != 0 else 50
            
            strand_info = {
                'StrandID': strand_id,
                'OligoKey': '.'.join([str(x) for x in oligo.key]),
                'HelixID': strand.vh,
                'StartIdx': strand.idx5p,
                'EndIdx': strand.idx3p,
                'Direction': 'fwd' if strand.direction > 0 else 'rev',
                'Length': len(sequence),
                'Sequence': sequence,
                'ProbFold': edge_prob,
                'LogProbFold': edge_logprob,
                'Tf': edge_Tf,
                'TfColor': 0,
                'maxTm': tm_value,
                'dGtotal': dG_total + dG_conc,  # Include concentration term, basically 3.225
                'dGhyb': dG_total,  # Hybridization free energy
                'dGloop': 0,  # No loop for individual strands
                'dGconc': dG_conc,
            }
            return strand_info
            
        except Exception as e:
            print(f"DEBUG: Error calculating thermodynamics for strand {strand_id}: {e}")
            import traceback
            traceback.print_exc()
    
    def _get_tf_color(self, tf, min_tf, max_tf):
        """Convert Tf value to hex color using coolwarm colormap."""
        import matplotlib.cm as cm
        
        cmap = cm.get_cmap('coolwarm').reversed()
        
        # Normalize Tf to [0, 1] range
        if tf <= min_tf:
            normalized_tf = 0
        elif tf >= max_tf:
            normalized_tf = 1
        else:
            normalized_tf = (tf - min_tf) / (max_tf - min_tf)
        
        rgb = cmap(normalized_tf)[:3]
        return to_hex(rgb)
    
    def create_strand_data_file(self):
        """Create Excel file with strand-level data."""
        strand_data = self.extract_strand_data()
        
        if not strand_data:
            print("No strand data found.")
            return None
        
        # Convert to DataFrame
        df = pd.DataFrame(strand_data)
        
        # Create Excel file
        with pd.ExcelWriter(self.strand_data_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='StrandData', index=False)
            
            # Add summary sheet
            summary_data = {
                'Total Strands': len(df),
                'Average Length': df['Length'].mean(),
                'Average Tf': df['Tf'].mean(),
                'Average ProbFold': df['ProbFold'].mean(),
            }
            summary_df = pd.DataFrame(list(summary_data.items()), columns=['Parameter', 'Value'])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        print(f"Created strand data file: {self.strand_data_file}")
        #NEW Moved dGloop calculation here
        from autobreak import utilities
        # Get scaffold info for dG_loop calculations
        scaffold = self.origami.oligos['scaffold'][0] if self.origami.oligos['scaffold'] else None
        scaffold_length = scaffold.length if scaffold else 1000
        is_circular = scaffold.circular if scaffold else True
        temperature_kelvin = self.autobreak.optim_temperature_kelvin
        
        # Group strands by oligo to find connections
        oligo_groups = {}
        for _, row in df.iterrows():
            oligo_key = row['OligoKey']
            if oligo_key not in oligo_groups:
                oligo_groups[oligo_key] = []
            oligo_groups[oligo_key].append(row)
        
        # Calculate dG_loop for each strand
        dg_loop_5p = {}  # 5' side (to previous strand)
        dg_loop_3p = {}  # 3' side (to next strand)
        
        for oligo_key, strands in oligo_groups.items():
            for i, strand in enumerate(strands):
                strand_id = strand['StrandID']
                
                # 3' dG_loop (to next strand)
                if i < len(strands) - 1:
                    next_strand = strands[i + 1]
                    current_strand_obj = self._find_strand_by_id_from_df(strand['StrandID'])
                    next_strand_obj = self._find_strand_by_id_from_df(next_strand['StrandID'])
                    
                    if current_strand_obj and next_strand_obj:
                        current_pos = self._get_mean_scaffold_pos(current_strand_obj)
                        next_pos = self._get_mean_scaffold_pos(next_strand_obj)
                        
                        if current_pos is not None and next_pos is not None:
                            dG_loop, _ = utilities.position_to_loop_dG(
                                scaffold_length, current_pos, next_pos,
                                is_circular, temperature_kelvin
                            )
                            dg_loop_3p[strand_id] = dG_loop
                        else:
                            dg_loop_3p[strand_id] = 0
                    else:
                        dg_loop_3p[strand_id] = 0
                else:
                    dg_loop_3p[strand_id] = 0
                
                # 5' dG_loop (from previous strand)
                if i > 0:
                    prev_strand = strands[i - 1]
                    prev_strand_obj = self._find_strand_by_id_from_df(prev_strand['StrandID'])
                    current_strand_obj = self._find_strand_by_id_from_df(strand['StrandID'])
                    
                    if prev_strand_obj and current_strand_obj:
                        prev_pos = self._get_mean_scaffold_pos(prev_strand_obj)
                        current_pos = self._get_mean_scaffold_pos(current_strand_obj)
                        
                        if prev_pos is not None and current_pos is not None:
                            dG_loop, _ = utilities.position_to_loop_dG(
                                scaffold_length, prev_pos, current_pos,
                                is_circular, temperature_kelvin
                            )
                            dg_loop_5p[strand_id] = dG_loop
                        else:
                            dg_loop_5p[strand_id] = 0
                    else:
                        dg_loop_5p[strand_id] = 0
                else:
                    dg_loop_5p[strand_id] = 0
        
        # Add dG_loop columns to dataframe
        df['dGloop_5p'] = df['StrandID'].map(dg_loop_5p).fillna(0)
        df['dGloop_3p'] = df['StrandID'].map(dg_loop_3p).fillna(0)
        return df
    
    def create_strand_svg_heatmap(self, df=None):
        """Generate strand-level SVG heatmap using cn2svg for base structure."""
        if df is None:
            df = self.create_strand_data_file()
        
        if df is None or df.empty:
            print("No data available for SVG creation.")
            return
        
        # Use cn2svg to generate the base structure (helix labels, grid, etc.)
        # then overlay strand-level rectangles
        success = self._generate_cn2svg_base()
        
        if success:
            # Overlay strand-level data on the cn2svg base
            self._overlay_strand_data_on_cn2svg(df)
        else:
            print("cn2svg generation failed")
        
        print(f"Created strand-level SVG heatmap: {self.strand_svg_path}")
    
    def _generate_cn2svg_base(self):
        """Generate the base SVG structure using cn2svg."""
        try:
            # Use cn2svg to create the structural elements (helix labels, grid, etc.)
            from cn2svg import cn2svg
            
            # Create document and generate base SVG
            if hasattr(self.autobreak, "json_legacy_output") and os.path.exists(self.autobreak.json_legacy_output):
                json_path = self.autobreak.json_legacy_output
            else:
                json_path = self.output_name
            cndoc = cn2svg.CadnanoDocument(json_path, "")

            # Generate path SVG with heatmap mode for structure
            path_svg_generator = cn2svg.CadnanoPathSvg(
                cndoc, 
                self.strand_svg_path,
                heatmap=True,  # This gives us the right structure
                scale=10
            )
            
            return True
            
        except Exception as e:
            print(f"Error generating cn2svg base: {e}")
            return False

    def _add_colored_crossovers_to_group(self, crossover_group, df, root):
        """Recreate crossover paths with individual colors based on dG_loop values."""
        import matplotlib.cm as cm
        from cn2svg import cn2svg
        
        # Collect all non-zero dG_loop values for normalization
        all_dg_values = []
        for _, row in df.iterrows():
            dg_3p = row['dGloop_3p']
            if dg_3p != 0:
                all_dg_values.append(dg_3p)
        
        if not all_dg_values:
            print("Failed to find any non-zero dG_loop values for crossovers.")
            return
        
        # Calculate min/max for normalization (only non-zero values)
        dg_min = min(all_dg_values)
        dg_max = max(all_dg_values)
        
        # Get colormap
        cmap = cm.get_cmap('coolwarm')
        crossover_group.clear()
        
        # Get layout parameters from cn2svg
        _BW = 10  # base_width (scale * vh_radius)
        _BH = 10  # base_height
        PATH_X_PADDING = 40
        path_radius_scaled = 36
        path_vh_margin = 30
        _pX = PATH_X_PADDING + path_radius_scaled * 3 + _BW/2 - 77.5
        
        # Create helix to y-coordinate mapping
        try:
            cndoc = cn2svg.CadnanoDocument(self.autobreak.json_legacy_output, '')
            helix_order = cndoc.vh_order
        except:
            helix_order = sorted(df['HelixID'].unique())
        
        id_coords = {}
        for i, helix_id in enumerate(helix_order):
            id_coords[helix_id] = PATH_X_PADDING + path_vh_margin * i
        
        # Group strands by oligo and sort
        oligo_groups = {}
        for _, row in df.iterrows():
            oligo_key = row['OligoKey']
            if oligo_key not in oligo_groups:
                oligo_groups[oligo_key] = []
            oligo_groups[oligo_key].append(row)
        
        for oligo_key in oligo_groups:
            oligo_groups[oligo_key].sort(key=lambda x: x['StrandID'])
        
        # Draw individual crossover curves
        crossover_count = 0
        
        for oligo_key, strands in oligo_groups.items():
            for i in range(len(strands) - 1):
                current_strand = strands[i]
                next_strand = strands[i + 1]
                
                # Check if this is actually a crossover (different helices)
                if current_strand['HelixID'] == next_strand['HelixID']:
                    continue
                
                # Get dG_loop value (should be same as next_strand's dGloop_5p)
                dg_loop = current_strand['dGloop_3p']
                
                if dg_loop == 0:
                    continue
                
                # Calculate positions
                curr_helix = current_strand['HelixID']
                next_helix = next_strand['HelixID']
                
                curr_idx3p = current_strand['EndIdx']
                next_idx5p = next_strand['StartIdx']
                
                curr_is_fwd = current_strand['Direction'] == 'fwd'
                next_is_fwd = next_strand['Direction'] == 'fwd'
                
                # Calculate end point of current strand
                prevX = _pX + curr_idx3p * _BW
                if curr_idx3p > current_strand['StartIdx']:  # forward
                    prevY = id_coords[curr_helix] - _BH/2
                else:  # reverse
                    prevY = id_coords[curr_helix] + _BH/2
                
                # Calculate start point of next strand
                x = _pX + next_idx5p * _BW
                if next_idx5p < next_strand['EndIdx']:  # forward
                    y = id_coords[next_helix] - _BH/2
                else:  # reverse
                    y = id_coords[next_helix] + _BH/2
                
                # Create quadratic bezier control point
                if next_is_fwd:
                    x1 = x + abs(y - prevY) * 0.03
                else:
                    x1 = x - abs(y - prevY) * 0.03
                y1 = (y + prevY) / 2
                
                # Create path with quadratic curve
                path_data = "M %s, %s Q %s %s, %s %s" % (prevX, prevY, x1, y1, x, y)
                
                # Normalize dG and get color
                if dg_max > dg_min:
                    normalized = (dg_loop - dg_min) / (dg_max - dg_min)
                else:
                    normalized = 0.5
                
                rgb = cmap(normalized)[:3]
                color_hex = '#{:02x}{:02x}{:02x}'.format(
                    int(rgb[0]*255),
                    int(rgb[1]*255),
                    int(rgb[2]*255)
                )
                
                # Create path element
                from pysvg.shape import Path as PysvgPath
                
                path_elem = ET.Element('path')
                path_elem.set('d', path_data)
                path_elem.set('fill', 'none')
                path_elem.set('stroke', color_hex)
                path_elem.set('stroke-width', '2')
                path_elem.set('stroke-linejoin', 'round')
                path_elem.set('id', f'crossover_{crossover_count}')
                path_elem.set('data-dg-loop', f'{dg_loop:.2f}')
                
                crossover_group.append(path_elem)
                crossover_count += 1
    
    def _add_colored_crossovers_and_endpoints(self, crossover_group, endpoints_group, df, root):
        """Recreate crossover paths and 3' endpoints with individual colors based on dG_loop values."""
        import matplotlib.cm as cm
        from cn2svg import cn2svg
        
        # Collect all non-zero dG_loop values for normalization
        all_dg_values = []
        for _, row in df.iterrows():
            dg_3p = row['dGloop_3p']
            dg_5p = row['dGloop_5p']
            if dg_3p != 0:
                all_dg_values.append(dg_3p)
            if dg_5p != 0:
                all_dg_values.append(dg_5p)
        
        if not all_dg_values:
            print("Failed to find any non-zero dG_loop values.")
            return
        
        # Calculate min/max for normalization (only non-zero values)
        dg_min = min(all_dg_values)
        dg_max = max(all_dg_values)
        
        # Get colormap
        cmap = cm.get_cmap('coolwarm')
        
        # Clear existing groups
        crossover_group.clear()
        endpoints_group.clear()
        
        # Get layout parameters from cn2svg
        _BW = 10  # base_width (scale * vh_radius)
        _BH = 10  # base_height
        PATH_X_PADDING = 40
        path_radius_scaled = 36
        path_vh_margin = 30
        _pX = PATH_X_PADDING + path_radius_scaled * 3 + _BW/2 - 77.5
        
        # Create helix to y-coordinate mapping
        try:
            cndoc = cn2svg.CadnanoDocument(self.autobreak.json_legacy_output, '')
            helix_order = cndoc.vh_order
        except:
            helix_order = sorted(df['HelixID'].unique())
        
        id_coords = {}
        for i, helix_id in enumerate(helix_order):
            id_coords[helix_id] = PATH_X_PADDING + path_vh_margin * i
        
        # Group strands by oligo and sort
        oligo_groups = {}
        for _, row in df.iterrows():
            oligo_key = row['OligoKey']
            if oligo_key not in oligo_groups:
                oligo_groups[oligo_key] = []
            oligo_groups[oligo_key].append(row)
        
        for oligo_key in oligo_groups:
            oligo_groups[oligo_key].sort(key=lambda x: x['StrandID'])
        
        # Draw individual endpoints first, then crossover curves
        crossover_count = 0
        endpoint_count = 0
        
        # First pass: draw all 3' endpoints
        for oligo_key, strands in oligo_groups.items():
            for i in range(len(strands)):
                current_strand = strands[i]
                
                # Handle 3' endpoint (if no dGloop_3p, this strand ends here)
                if current_strand['dGloop_3p'] == 0:
                    color_hex = "#7E7E7E"  # Grey for isolated strands
                    
                    # Calculate 3' endpoint position
                    helix_id = current_strand['HelixID']
                    idx3p = current_strand['EndIdx']
                    is_fwd = current_strand['Direction'] == 'fwd'
                    
                    x = _pX + idx3p * _BW
                    
                    if is_fwd:
                        y = id_coords[helix_id] - _BH
                        x1, y1 = x, y
                        x2, y2 = x1 + _BW * 0.75, y1 + _BH / 2
                        x3, y3 = x, y1 + _BH
                    else:
                        y = id_coords[helix_id]
                        x1, y1 = x, y
                        x2, y2 = x - _BW * 0.75, y1 + _BH / 2
                        x3, y3 = x, y + _BH
                    
                    # Create polygon for 3' endpoint
                    pts = "%s,%s %s,%s %s,%s %s,%s" % (x1, y1, x2, y2, x3, y3, x1, y1)
                    polygon_elem = ET.Element('polygon')
                    polygon_elem.set('points', pts)
                    polygon_elem.set('fill', color_hex)
                    polygon_elem.set('stroke', 'none')
                    polygon_elem.set('stroke-width', '0.5')
                    polygon_elem.set('id', f'end3p_{endpoint_count}')
                    
                    endpoints_group.append(polygon_elem)
                    endpoint_count += 1
        
        # Second pass: draw all crossovers
        for oligo_key, strands in oligo_groups.items():
            for i in range(len(strands)):
                current_strand = strands[i]
                
                # Handle crossover (3' connection to next strand)
                if i < len(strands) - 1:
                    next_strand = strands[i + 1]
                    
                    # Check if this is actually a crossover (different helices)
                    if current_strand['HelixID'] != next_strand['HelixID']:
                        # Get dG_loop value
                        dg_loop = current_strand['dGloop_3p']
                        
                        if dg_loop != 0:
                            # Calculate positions
                            curr_helix = current_strand['HelixID']
                            next_helix = next_strand['HelixID']
                            
                            curr_idx3p = current_strand['EndIdx']
                            next_idx5p = next_strand['StartIdx']
                            
                            next_is_fwd = next_strand['Direction'] == 'fwd'
                            
                            # Calculate end point of current strand
                            prevX = _pX + curr_idx3p * _BW
                            if curr_idx3p > current_strand['StartIdx']:  # forward
                                prevY = id_coords[curr_helix] - _BH/2
                            else:  # reverse
                                prevY = id_coords[curr_helix] + _BH/2
                            
                            # Calculate start point of next strand
                            x = _pX + next_idx5p * _BW
                            if next_idx5p < next_strand['EndIdx']:  # forward
                                y = id_coords[next_helix] - _BH/2
                            else:  # reverse
                                y = id_coords[next_helix] + _BH/2
                            
                            # Create quadratic bezier control point
                            if next_is_fwd:
                                x1 = x + abs(y - prevY) * 0.03
                            else:
                                x1 = x - abs(y - prevY) * 0.03
                            y1 = (y + prevY) / 2
                            
                            # Create path with quadratic curve
                            path_data = "M %s, %s Q %s %s, %s %s" % (prevX, prevY, x1, y1, x, y)
                            
                            # Normalize dG and get color
                            if dg_max > dg_min:
                                normalized = (dg_loop - dg_min) / (dg_max - dg_min)
                            else:
                                normalized = 0.5
                            
                            rgb = cmap(normalized)[:3]
                            color_hex = '#{:02x}{:02x}{:02x}'.format(
                                int(rgb[0]*255),
                                int(rgb[1]*255),
                                int(rgb[2]*255)
                            )
                            
                            # Create path element
                            path_elem = ET.Element('path')
                            path_elem.set('d', path_data)
                            path_elem.set('fill', 'none')
                            path_elem.set('stroke', color_hex)
                            path_elem.set('stroke-width', '2')
                            path_elem.set('stroke-linejoin', 'round')
                            path_elem.set('id', f'crossover_{crossover_count}')
                            path_elem.set('data-dg-loop', f'{dg_loop:.2f}')
                            
                            crossover_group.append(path_elem)
                            crossover_count += 1
        
        print(f"Created {crossover_count} colored crossovers and {endpoint_count} colored 3' endpoints")

    def _overlay_strand_data_on_cn2svg(self, df):
        """Overlay individual strand rectangles on the cn2svg-generated base."""
        try:
            # Parse the cn2svg-generated SVG
            SVG_NS = "http://www.w3.org/2000/svg"
            ET.register_namespace("", SVG_NS)
            
            tree = ET.parse(self.strand_svg_path)
            root = tree.getroot()
            
            # Remove the original oligo paths since we want strand-level granularity
            oligos_group = None
            otherElements = None
            for element in root.iter():
                if element.get('id') == 'Oligos':
                    oligos_group = element
                    break
                elif element.get('id') == 'Endpoints':
                    otherElements = element
                    break
            
            if otherElements is not None:
                otherElements.clear()
                
            if oligos_group is not None:
                # Clear the oligo paths
                oligos_group.clear()
                oligos_group.attrib['id'] = 'StrandSegments'
                
                # Add strand-level rectangles
                self._add_strand_rectangles_to_group(oligos_group, df)
            
            # Find and modify crossover group in the same tree
            # crossover_group = None
            # for element in root.iter():
            #     if element.get('id') == 'OligoCrossovers':
            #         crossover_group = element
            #         break
            
            # if crossover_group is not None:
            #     # Add colored crossovers to the group
            #     self._add_colored_crossovers_to_group(crossover_group, df, root)
            # Find both crossover and endpoints groups
            crossover_group = None
            endpoints_group = None
            for element in root.iter():
                if element.get('id') == 'OligoCrossovers':
                    crossover_group = element
                elif element.get('id') == 'Endpoints':
                    endpoints_group = element
            
            if crossover_group is not None and endpoints_group is not None:
                for parent in root.iter():
                    if crossover_group in list(parent) and endpoints_group in list(parent):
                        # Get indices
                        children = list(parent)
                        crossover_idx = children.index(crossover_group)
                        endpoint_idx = children.index(endpoints_group)
                        
                        # Remove and re-insert to change order
                        parent.remove(endpoints_group)
                        parent.remove(crossover_group)
                        parent.insert(crossover_idx, endpoints_group)
                        parent.insert(endpoint_idx, crossover_group)
                        break
                # Add colored crossovers and endpoints
                self._add_colored_crossovers_and_endpoints(crossover_group, endpoints_group, df, root)
            
            
            # Add interactivity CSS
            style_element = ET.Element(f"{{{SVG_NS}}}style")
            style_element.text = """
            .strand-item {
                cursor: pointer;
                transition: all 0.2s ease;
            }
            .strand-item:hover {
                stroke: #000000;
                stroke-width: 2;
                opacity: 0.8;
            }
            .strand-item.highlight {
                stroke: #000000;
                stroke-width: 3;
                stroke-dasharray: 5,5;
                filter: brightness(1.2);
            }
            """
            root.insert(0, style_element)
            
            # Save the enhanced SVG
            tree.write(self.strand_svg_path, encoding='utf-8', xml_declaration=True)
            
            print(f"Overlaid {len(df)} strand segments on cn2svg base")
            
        except Exception as e:
            print(f"Error overlaying strand data: {e}")
            import traceback
            traceback.print_exc()
    
    def _add_strand_rectangles_to_group(self, parent_group, df):
        """Add individual strand rectangles to the SVG group with 7-base grid."""
        # Calculate positioning based on cn2svg's path view layout
        PATH_X_PADDING = 40
        path_radius_scaled = 36  # Approximate from cn2svg
        base_width = 10  # cn2svg scale
        base_height = 10
        
        x_offset = PATH_X_PADDING + path_radius_scaled
        
        # Group strands by helix to match cn2svg layout
        helix_groups = df.groupby('HelixID')
        
        # Get the helix order from the original design
        try:
            from cn2svg import cn2svg
            #cndoc = cn2svg.CadnanoDocument(self.autobreak.json_legacy_output, '')
            if os.path.exists(self.autobreak.json_legacy_output):
                json_path = self.autobreak.json_legacy_output
            else:
                json_path = self.output_name
            cndoc = cn2svg.CadnanoDocument(json_path, "")
            helix_order = cndoc.vh_order
        except:
            helix_order = sorted(helix_groups.groups.keys())
        
        path_vh_margin = 30  # Spacing between helices in cn2svg
        
        # First, add grid lines with 7-base intervals
        min_pos = df[['StartIdx', 'EndIdx']].values.min()
        max_pos = df[['StartIdx', 'EndIdx']].values.max()
        
        grid_group = ET.Element('g')
        grid_group.set('id', 'grid-lines')
        grid_group.set('stroke', '#cccccc')
        grid_group.set('stroke-width', '0.5')
        grid_group.set('opacity', '0.7')
        
        # Create vertical grid lines every 7 bases
        for pos in range(int(min_pos), int(max_pos) + 1, 7):
            x_coord = x_offset + pos * base_width - 15 
            
            # Draw vertical line spanning all helices
            line = ET.Element('line')
            line.set('x1', str(x_coord))
            line.set('y1', str(PATH_X_PADDING))
            line.set('x2', str(x_coord))
            line.set('y2', str(PATH_X_PADDING + len(helix_order) * path_vh_margin))
            line.set('stroke', "#838383")
            grid_group.append(line)
        
        parent_group.append(grid_group)

        # Add scaffold lines by tracing scaffold oligo v2
        # Add scaffold lines - continuous path through all strands
        scaffold_group = ET.Element('g')
        scaffold_group.set('id', 'scaffold-lines')
        scaffold_group.set('stroke', '#0066cc')
        scaffold_group.set('stroke-width', '2')
        scaffold_group.set('fill', 'none')

        # Get scaffold oligo
        scaffold_oligo = self.origami.oligos['scaffold'][0] if self.origami.oligos['scaffold'] else None

        if scaffold_oligo:
            current_strand = scaffold_oligo.null_strand.next_strand
            path_segments = []
            
            while current_strand:
                vh = current_strand.vh
                
                if vh in helix_order:
                    helix_idx = helix_order.index(vh)
                    if helix_idx % 2 == 0:
                        y_center = PATH_X_PADDING + path_vh_margin * helix_idx - 5
                    else:
                        y_center = PATH_X_PADDING + path_vh_margin * helix_idx + 5
                    
                    idx5p = current_strand.idx5p
                    idx3p = current_strand.idx3p
                    
                    # Scaffold follows strand direction
                    x_start = x_offset + idx5p * base_width
                    x_end = x_offset + idx3p * base_width
                    
                    if not path_segments:
                        path_segments.append(f"M {x_start} {y_center}")
                    else:
                        path_segments.append(f"L {x_start} {y_center}")
                    
                    path_segments.append(f"L {x_end} {y_center}")
                
                current_strand = current_strand.next_strand
            
            if path_segments:
                scaffold_path = ET.Element('path')
                scaffold_path.set('d', ' '.join(path_segments))
                scaffold_group.append(scaffold_path)

        parent_group.append(scaffold_group)

        # Initialize strand connections list
        strand_connections = []
        
        # Now add strand rectangles
        for helix_idx, helix_id in enumerate(helix_order):
            if helix_id not in helix_groups.groups:
                continue
                
            helix_strands = helix_groups.get_group(helix_id)
            y_base = PATH_X_PADDING + path_vh_margin * helix_idx
            
            for _, strand in helix_strands.iterrows():
                # Calculate strand position
                start_pos = min(strand['StartIdx'], strand['EndIdx'])
                end_pos = max(strand['StartIdx'], strand['EndIdx'])
                length = end_pos - start_pos
                
                x = x_offset + start_pos * base_width
                width = length * base_width
                
                # Position based on direction (forward = top, reverse = bottom)
                if strand['Direction'] == 'fwd':
                    y = y_base - base_height
                    height = base_height
                else:
                    y = y_base
                    height = base_height
                
                # Create rectangle element
                rect_elem = ET.Element('rect')
                rect_elem.set('x', str(x))
                rect_elem.set('y', str(y))
                rect_elem.set('width', str(width))
                rect_elem.set('height', str(height))
                rect_elem.set('fill', strand['TfColor'])
                rect_elem.set('stroke', strand['TfColor'])
                rect_elem.set('stroke-width', '0.5')
                rect_elem.set('class', 'strand-item')
                rect_elem.set('data-strand-id', str(strand['StrandID']))
                rect_elem.set('data-helix-id', str(strand['HelixID']))
                rect_elem.set('data-tf', str(strand['Tf']))
                rect_elem.set('data-length', str(strand['Length']))
                rect_elem.set('data-sequence', strand['Sequence'])
                
                parent_group.append(rect_elem)

                # Store connection info
                strand_connections.append({
                    'strand_id': strand['StrandID'],
                    'oligo_key': strand['OligoKey'],
                    'start_x': x,
                    'end_x': x + width,
                    'y': y_base,
                    'direction': strand['Direction']
                })

    
    def create_interactive_html(self, df=None):
        """Create interactive HTML combining SVG heatmap with clickable strand data."""
        if df is None:
            df = pd.read_excel(self.strand_data_file, sheet_name='StrandData') if os.path.exists(self.strand_data_file) else None
        
        if df is None or df.empty:
            print("No strand data available for interactive HTML.")
            return
        
        # Ensure SVG exists
        if not os.path.exists(self.strand_svg_path):
            self.create_strand_svg_heatmap(df)
        
        # Process SVG to add interactivity
        svg_content = self._process_svg_for_interactivity(df)
        
        # Create HTML table
        table_html = self._create_data_table(df)
        
        # Generate complete HTML
        html_content = self._generate_html_template(svg_content, table_html)
        
        # Write HTML file
        with open(self.interactive_html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"Created interactive HTML: {self.interactive_html_path}")
    
    def _process_svg_for_interactivity(self, df):
        """Process SVG to add click handlers and IDs."""
        if not os.path.exists(self.strand_svg_path):
            return ""
        
        # Parse SVG
        SVG_NS = "http://www.w3.org/2000/svg"
        ET.register_namespace("", SVG_NS)
        
        tree = ET.parse(self.strand_svg_path)
        root = tree.getroot()
        
        # Find rectangles and add interactivity based on their gid attributes
        strand_count = 0
        for element in root.iter():
            if 'rect' in element.tag.lower():
                gid = element.get('gid')
                if gid and gid.startswith('strand_'):
                    # This rectangle represents a strand
                    element.set("class", "strand-item")
                    element.set("data-key", gid.replace('strand_', ''))
                    element.set("style", "cursor: pointer; transition: stroke-width 0.2s ease;")
                    strand_count += 1
        
        # Add enhanced CSS for interactivity
        style_element = ET.Element(f"{{{SVG_NS}}}style")
        style_element.text = """
        .strand-item {
            cursor: pointer !important;
            transition: all 0.2s ease;
        }
        .strand-item:hover {
            stroke: #333333 !important;
            stroke-width: 2 !important;
            opacity: 0.8;
        }
        .strand-item.highlight {
            stroke: #000000 !important;
            stroke-width: 3 !important;
            stroke-dasharray: 5,5 !important;
        }
        text {
            pointer-events: none;
            user-select: none;
        }
        circle {
            pointer-events: none;
        }
        line {
            pointer-events: none;
        }
        """
        root.insert(0, style_element)
        return ET.tostring(root, encoding="unicode", method="xml")
    
    def _create_data_table(self, df):
        """Create HTML table with sortable columns and heatmap coloring."""

        display_cols = ['StrandID', 'OligoKey', 'HelixID', 'StartIdx', 'EndIdx', 
                    'Direction', 'Length', 'Sequence', 'Tf', 'ProbFold', 'dGtotal',
                    'dGloop_5p', 'dGloop_3p']
        
        available_cols = [col for col in display_cols if col in df.columns]
        df_display = df[available_cols].copy()
        
        numeric_cols = df_display.select_dtypes(include=[np.number]).columns
        df_display[numeric_cols] = df_display[numeric_cols].round(3)
        
        # Create table HTML with sortable headers
        thead = "<tr>" + "".join(
            f'<th class="sortable" onclick="sortTable({i})">{html.escape(str(col))} <span class="sort-arrow">↕</span></th>' 
            for i, col in enumerate(df_display.columns)
        ) + "</tr>"
        
        # Get min/max for heatmap coloring
        tf_min, tf_max = df_display['Tf'].min(), df_display['Tf'].max()
        
        rows = []
        for i, row in df_display.reset_index(drop=True).iterrows():
            tds = []
            for col in df_display.columns:
                value = row[col]
                
                # Apply heatmap coloring ONLY to Tf
                bg_color = ''
                if col == 'Tf':
                    normalized = (value - tf_min) / (tf_max - tf_min) if tf_max > tf_min else 0
                    import matplotlib.cm as cm
                    cmap = cm.get_cmap('coolwarm').reversed()
                    rgb = cmap(normalized)[:3]
                    bg_color = f'background-color: rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)});'
                
                tds.append(f'<td style="{bg_color}">{html.escape(str(value))}</td>')
            
            rows.append(f'<tr class="row" data-key="{i}">{"".join(tds)}</tr>')
        
        return f'<table id="strand-table"><thead>{thead}</thead><tbody>{"".join(rows)}</tbody></table>'
    
    def _find_strand_by_id_from_df(self, strand_id):
        """Find actual strand object by ID."""
        current_id = 0
        for oligo in self.origami.oligos['staple']:
            current_strand = oligo.null_strand.next_strand
            while current_strand:
                if current_id == strand_id:
                    return current_strand
                current_id += 1
                current_strand = current_strand.next_strand
        return None

    def _get_mean_scaffold_pos(self, strand):
        """Get mean scaffold position for a strand."""
        if not hasattr(strand, 'scaffoldPos') or not strand.scaffoldPos:
            return None
        valid_pos = [p for p in strand.scaffoldPos if p is not None]
        if not valid_pos:
            return None
        return sum(valid_pos) / len(valid_pos)

    def _generate_html_template(self, svg_content, table_html):
        """Generate complete HTML template with proper syntax."""
        # Enhanced CSS with persistent grid labels
        css = """
        body { 
            margin: 0; 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
            overflow: hidden;
        }
        #app { 
            display: grid; 
            grid-template-columns: 2fr 1fr; 
            height: 100vh; 
        }
        #left, #right { 
            overflow: hidden; 
            display: flex;
            flex-direction: column;
        }
        #left { 
            border-right: 2px solid #ddd; 
            background: #f8f9fa; 
        }
        #svg-container { 
            flex: 1;
            overflow: auto;
            position: relative;
            background: white;
        }
        #svg-container svg { 
            display: block; 
            cursor: crosshair;
            transition: transform 0.2s ease;
        }
        .grid-labels {
            position: fixed;
            top: 80px;
            left: 0;
            right: 50%;
            height: 30px;
            background: rgba(248, 249, 250, 0.95);
            border-bottom: 1px solid #ddd;
            z-index: 100;
            pointer-events: none;
            display: flex;
            align-items: center;
            font-size: 10px;
            color: #666;
            overflow: hidden;
        }
        .grid-label {
            position: absolute;
            transform: translateX(-50%);
        }
        .sortable {
            cursor: pointer;
            user-select: none;
        }
        .sortable:hover {
            background: #e0e0e0;
        }
        .sort-arrow {
            font-size: 10px;
            color: #999;
            margin-left: 4px;
        }
        #strand-table-container {
            flex: 1;
            overflow: auto;
        }
        #strand-table { 
            width: 100%; 
            border-collapse: collapse; 
            font-size: 11px; 
        }
        #strand-table th, #strand-table td { 
            border-bottom: 1px solid #eee; 
            padding: 4px 6px; 
            text-align: left; 
        }
        #strand-table th { 
            background: #f5f5f5; 
            font-weight: 600; 
            position: sticky; 
            top: 0; 
            z-index: 5;
        }
        #strand-table tr:hover { 
            background: #f0f8ff; 
        }
        #strand-table tr.selected { 
            background: #fff3cd !important; 
            outline: 2px solid #ffc107;
        }
        .controls { 
            padding: 10px; 
            border-bottom: 2px solid #eee; 
            background: #f8f9fa; 
            position: sticky; 
            top: 0; 
            z-index: 10; 
        }
        .controls h3 { 
            margin: 0 0 8px 0; 
            color: #333; 
        }
        .controls-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 15px;
            font-size: 11px;
            color: #666;
        }
        .zoom-controls {
            display: flex;
            gap: 5px;
            align-items: center;
        }
        .zoom-btn {
            background: #007bff;
            color: white;
            border: none;
            padding: 5px 10px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 12px;
        }
        .zoom-btn:hover {
            background: #0056b3;
        }
        .zoom-btn:disabled {
            background: #6c757d;
            cursor: not-allowed;
        }
        .zoom-info {
            font-size: 11px;
            color: #666;
            margin-left: 10px;
        }
        .nav-btn {
        background: #007bff;
        color: white;
        border: none;
        padding: 2px 8px;
        border-radius: 3px;
        cursor: pointer;
        font-size: 14px;
        }
        .nav-btn:hover {
            background: #0056b3;
        }
        .legend {
            font-size: 11px;
            color: #666;
            margin-top: 5px;
        }
        """
        
        js = """
        let currentZoom = 1;
        let currentTransform = { x: 0, y: 0 };
        let selectedStrandId = null;
        let isDragging = false;
        let dragStart = { x: 0, y: 0 };
        let svgElement = null;
        let svgContainer = null;
        let gridLabelsContainer = null;
        
        function navigateToStrand(strandId) {
            clearHighlights();
            
            // Highlight table row
            const row = document.querySelector(`tr[data-key="${strandId}"]`);
            if (row) {
                row.classList.add('selected');
                row.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
            
            // Highlight and center SVG element
            const strandElement = svgElement.querySelector(`rect[data-strand-id="${strandId}"], .strand-item[data-strand-id="${strandId}"]`);
            if (strandElement) {
                strandElement.classList.add('highlight');
                strandElement.style.stroke = '#000000';
                strandElement.style.strokeWidth = '3';
                strandElement.style.strokeDasharray = '5,5';
                
                // Get strand position in SVG coordinates
                const strandRect = strandElement.getBBox();
                const strandCenterX = strandRect.x + strandRect.width / 2;
                const strandCenterY = strandRect.y + strandRect.height / 2;
                
                // Get container dimensions (left half of screen)
                const containerRect = svgContainer.getBoundingClientRect();
                const containerCenterX = containerRect.width / 2;
                const containerCenterY = containerRect.height / 2;
                
                // Calculate new transform to center the strand
                currentTransform.x = containerCenterX - strandCenterX * currentZoom;
                currentTransform.y = containerCenterY - strandCenterY * currentZoom;
                
                // Apply transform with smooth transition
                svgElement.style.transition = 'transform 0.5s ease';
                updateTransform();
                
                // Remove transition after animation
                setTimeout(() => {
                    svgElement.style.transition = '';
                }, 500);
            }
            
            selectedStrandId = strandId;
        }

        let sortDirections = {};

        function sortTable(columnIndex) {
            const table = document.getElementById('strand-table');
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr.row'));
            
            // Toggle sort direction
            if (!sortDirections[columnIndex]) {
                sortDirections[columnIndex] = 'asc';
            } else {
                sortDirections[columnIndex] = sortDirections[columnIndex] === 'asc' ? 'desc' : 'asc';
            }
            
            const direction = sortDirections[columnIndex];
            const headers = Array.from(table.querySelectorAll('th.sortable'));
            const columnName = headers[columnIndex].textContent.trim().replace('↕', '').replace('↑', '').replace('↓', '').trim();
            
            rows.sort((a, b) => {
                let comparison = 0;
                
                // Special handling for OligoKey column
                if (columnName === 'OligoKey') {
                    const aOligo = a.cells[columnIndex].textContent.trim();
                    const bOligo = b.cells[columnIndex].textContent.trim();
                    const aHelix = parseInt(a.cells[columnIndex + 1].textContent.trim());
                    const bHelix = parseInt(b.cells[columnIndex + 1].textContent.trim());
                    
                    // Parse OligoKey as "num1.num2.-1"
                    const aParts = aOligo.split('.');
                    const bParts = bOligo.split('.');
                    
                    const aFirst = parseInt(aParts[0]);
                    const bFirst = parseInt(bParts[0]);
                    const aSecond = parseInt(aParts[1]);
                    const bSecond = parseInt(bParts[1]);
                    
                    // First sort by first number
                    if (aFirst !== bFirst) {
                        comparison = aFirst - bFirst;
                    }
                    // Then by second number
                    else if (aSecond !== bSecond) {
                        comparison = aSecond - bSecond;
                    }
                    // Finally by HelixID
                    else {
                        comparison = aHelix - bHelix;
                    }
                } else {
                    const aText = a.cells[columnIndex].textContent.trim();
                    const bText = b.cells[columnIndex].textContent.trim();
                    
                    // Try to parse as numbers
                    const aNum = parseFloat(aText);
                    const bNum = parseFloat(bText);
                    
                    if (!isNaN(aNum) && !isNaN(bNum)) {
                        comparison = aNum - bNum;
                    } else {
                        comparison = aText.localeCompare(bText);
                    }
                }
                
                return direction === 'asc' ? comparison : -comparison;
            });
            
            // Re-append sorted rows
            rows.forEach(row => tbody.appendChild(row));
            
            // Update sort arrows
            headers.forEach((header, i) => {
                const arrow = header.querySelector('.sort-arrow');
                if (i === columnIndex) {
                    arrow.textContent = direction === 'asc' ? '↑' : '↓';
                } else {
                    arrow.textContent = '↕';
                }
            });
        }

        function initializeViewer() {
            svgElement = document.querySelector('#svg-container svg');
            svgContainer = document.getElementById('svg-container');
            
            // Create persistent grid labels container
            gridLabelsContainer = document.createElement('div');
            gridLabelsContainer.className = 'grid-labels';
            gridLabelsContainer.id = 'grid-labels';
            document.body.appendChild(gridLabelsContainer);
            
            if (svgElement) {
                // Add wheel zoom
                svgContainer.addEventListener('wheel', handleWheel, { passive: false });
                
                // Add mouse drag for panning
                svgContainer.addEventListener('mousedown', handleMouseDown);
                svgContainer.addEventListener('mousemove', handleMouseMove);
                svgContainer.addEventListener('mouseup', handleMouseUp);
                svgContainer.addEventListener('mouseleave', handleMouseUp);
                
                // Prevent default drag behavior
                svgContainer.addEventListener('dragstart', e => e.preventDefault());
                
                // Initialize strand click handlers
                initializeStrandClickHandlers();
                
                // Update grid labels initially
                updateGridLabels();
                updateZoomInfo();
            }
        }
        
        function updateGridLabels() {
            if (!gridLabelsContainer || !svgContainer) return;
            
            // Clear existing labels
            gridLabelsContainer.innerHTML = '';
            
            // Get container dimensions
            const containerRect = svgContainer.getBoundingClientRect();
            const containerWidth = containerRect.width;
            
            // Calculate visible range based on current transform and zoom
            const visibleStartX = -currentTransform.x / currentZoom;
            const visibleEndX = (containerWidth - currentTransform.x) / currentZoom;
            
            // Grid interval every 7 bases (matching SVG)
            const gridInterval = 7;
            const baseWidth = 10; // Match cn2svg base width
            
            // Calculate start and end positions in base units
            // Constants matching SVG generation  
            const PATH_X_PADDING = 40;
            const path_radius_scaled = 36;
            const x_offset = PATH_X_PADDING + path_radius_scaled;

            // Calculate start and end positions accounting for x_offset
            const startPos = Math.floor((visibleStartX - x_offset) / baseWidth / gridInterval) * gridInterval;
            const endPos = Math.ceil((visibleEndX - x_offset) / baseWidth / gridInterval) * gridInterval;
            // Create labels for visible grid positions
            for (let pos = startPos; pos <= endPos; pos += gridInterval) {
                if (pos < 0) continue;
                
                const pixelX = (x_offset + pos * baseWidth) * currentZoom + currentTransform.x;

                // Only show labels that are within the visible area
                if (pixelX >= -50 && pixelX <= containerWidth + 50) {
                    const label = document.createElement('div');
                    label.className = 'grid-label';
                    label.textContent = pos.toString();
                    label.style.left = pixelX + 'px';
                    gridLabelsContainer.appendChild(label);
                }
            }
        }
        
        function initializeStrandClickHandlers() {
            // Find all strand elements
            const strandElements = svgElement.querySelectorAll('.strand-item, rect[data-strand-id]');
            
            strandElements.forEach((element) => {
                const strandId = element.getAttribute('data-strand-id');
                
                if (strandId) {
                    // Add click handler
                    element.addEventListener('click', function(e) {
                        e.stopPropagation();
                        highlightStrand(this, strandId);
                    });
                    
                    // Add hover effect
                    element.addEventListener('mouseenter', function() {
                        if (!this.classList.contains('highlight')) {
                            this.style.stroke = '#333333';
                            this.style.strokeWidth = '2';
                        }
                    });
                    
                    element.addEventListener('mouseleave', function() {
                        if (!this.classList.contains('highlight')) {
                            this.style.stroke = '#000000';
                            this.style.strokeWidth = '0.5';
                        }
                    });
                }
            });
            
            console.log(`Initialized ${strandElements.length} strand click handlers`);
        }
        
        function highlightStrand(element, strandId) {
            // Clear previous highlights
            clearHighlights();
            
            // Highlight clicked strand
            element.classList.add('highlight');
            element.style.stroke = '#000000';
            element.style.strokeWidth = '3';
            element.style.strokeDasharray = '5,5';
            selectedStrandId = strandId;
            
            // Highlight corresponding table row
            const row = document.querySelector(`tr[data-key="${strandId}"]`);
            if (row) {
                row.classList.add('selected');
                row.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
            
            console.log(`Highlighted strand ${strandId}`);
        }
        
        function clearHighlights() {
            // Clear SVG highlights
            document.querySelectorAll('.strand-item.highlight, rect.highlight').forEach(el => {
                el.classList.remove('highlight');
                el.style.stroke = '#000000';
                el.style.strokeWidth = '0.5';
                el.style.strokeDasharray = '';
            });
            
            // Clear table highlights
            document.querySelectorAll('tr.selected').forEach(tr => tr.classList.remove('selected'));
        }
        
        function handleWheel(e) {
            e.preventDefault();
            
            const rect = svgContainer.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            
            const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
            const newZoom = Math.min(Math.max(currentZoom * zoomFactor, 0.1), 10);
            
            if (newZoom !== currentZoom) {
                const zoomPointX = (mouseX - currentTransform.x) / currentZoom;
                const zoomPointY = (mouseY - currentTransform.y) / currentZoom;
                
                currentTransform.x = mouseX - zoomPointX * newZoom;
                currentTransform.y = mouseY - zoomPointY * newZoom;
                currentZoom = newZoom;
                
                updateTransform();
                updateGridLabels();
            }
        }
        
        function handleMouseDown(e) {
            if (e.button === 0) {
                isDragging = true;
                dragStart.x = e.clientX - currentTransform.x;
                dragStart.y = e.clientY - currentTransform.y;
                svgContainer.style.cursor = 'grabbing';
            }
        }
        
        function handleMouseMove(e) {
            if (isDragging) {
                currentTransform.x = e.clientX - dragStart.x;
                currentTransform.y = e.clientY - dragStart.y;
                updateTransform();
                updateGridLabels();
            }
        }
        
        function handleMouseUp() {
            isDragging = false;
            svgContainer.style.cursor = 'crosshair';
        }
        
        function updateTransform() {
            if (svgElement) {
                svgElement.style.transform = `translate(${currentTransform.x}px, ${currentTransform.y}px) scale(${currentZoom})`;
                svgElement.style.transformOrigin = '0 0';
                updateZoomInfo();
            }
        }
        
        function zoomIn() {
            const newZoom = Math.min(currentZoom * 1.2, 10);
            if (newZoom !== currentZoom) {
                currentZoom = newZoom;
                updateTransform();
                updateGridLabels();
            }
        }
        
        function zoomOut() {
            const newZoom = Math.max(currentZoom * 0.8, 0.1);
            if (newZoom !== currentZoom) {
                currentZoom = newZoom;
                updateTransform();
                updateGridLabels();
            }
        }
        
        function resetZoom() {
            currentZoom = 1;
            currentTransform = { x: 0, y: 0 };
            updateTransform();
            updateGridLabels();
        }
        
        function updateZoomInfo() {
            const zoomInfo = document.querySelector('.zoom-info');
            if (zoomInfo) {
                zoomInfo.textContent = `Zoom: ${(currentZoom * 100).toFixed(0)}%`;
            }
            
            const zoomOutBtn = document.getElementById('zoom-out');
            const zoomInBtn = document.getElementById('zoom-in');
            if (zoomOutBtn) zoomOutBtn.disabled = currentZoom <= 0.1;
            if (zoomInBtn) zoomInBtn.disabled = currentZoom >= 10;
        }
        
        // Table row click handler
        document.addEventListener('DOMContentLoaded', function() {
            initializeViewer();
            
            const tbody = document.querySelector('#strand-table tbody');
            if (tbody) {
                tbody.addEventListener('click', function(e) {
                    const tr = e.target.closest('tr.row');
                    if (!tr) return;
                    
                    const strandId = tr.getAttribute('data-key');
                    clearHighlights();
                    
                    // Highlight table row
                    tr.classList.add('selected');
                    
                    // Find and highlight corresponding strand in SVG
                    const strandElement = svgElement.querySelector(`rect[data-strand-id="${strandId}"], .strand-item[data-strand-id="${strandId}"]`);
                    
                    if (strandElement) {
                        strandElement.classList.add('highlight');
                        strandElement.style.stroke = '#000000';
                        strandElement.style.strokeWidth = '3';
                        strandElement.style.strokeDasharray = '5,5';
                        // strandElement.scrollIntoView({behavior: 'smooth', block: 'center'});
                        navigateToStrand(strandId);
                        updateTransform();
                        updateGridLabels();
                    }
                    
                    selectedStrandId = parseInt(strandId);
                    console.log(`Selected strand ${strandId} from table`);
                });
            }
        });
        """
        
        # Create the final HTML template
        html_template = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Interactive Strand-Level DNA Origami Visualization</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>{}</style>
</head>
<body>
    <div id="app">
        <div id="left">
            <div class="controls">
                <h3>Interactive Strand Heatmap</h3>
                <div class="controls-grid">
                    <div>
                        <div class="zoom-controls">
                            <button id="zoom-out" class="zoom-btn" onclick="zoomOut()">-</button>
                            <button id="zoom-in" class="zoom-btn" onclick="zoomIn()">+</button>
                            <button class="zoom-btn" onclick="resetZoom()">Reset</button>
                            <span class="zoom-info">Zoom: 100%</span>
                        </div>
                        <div class="legend">Drag to pan • Scroll to zoom • Click strands for details</div>
                    </div>
                    <div>
                        <div>Individual strand segments with thermodynamic data</div>
                        <div>Colors based on folding temperature (Tf)</div>
                    </div>
                </div>
            </div>
            <div id="svg-container">{}</div>
        </div>
        <div id="right">
            <div class="controls">
                <h3>Strand Properties</h3>
                <div style="font-size: 11px; color: #666;">
                    Click table rows or strands to highlight
                </div>
            </div>
            <div id="strand-table-container">{}</div>
        </div>
    </div>
    <script>{}</script>
</body>
</html>""".format(css, svg_content, table_html, js)
        
        return html_template
    
    def generate_all_visualizations(self):
        """Generate all strand-level visualizations."""
        print("Generating strand-level visualizations...")
        
        # 1. Extract strand data
        df = self.create_strand_data_file()
        
        # 2. Create SVG heatmap
        self.create_strand_svg_heatmap(df)
        
        # 3. Create interactive HTML
        self.create_interactive_html(df)
        
        print("Strand visualization generation complete!")
        
        return {
            'strand_data_file': self.strand_data_file,
            'strand_svg_path': self.strand_svg_path,
            'interactive_html_path': self.interactive_html_path
        }


def create_strand_visualizations(autobreak_instance, override):
    """
    Convenience function to create strand-level visualizations.
    
    Args:
        autobreak_instance: Completed AutoBreak instance
        Override: Whether this is an initial run (pre-autobreak) or after autobreak processing.
        
    Returns:
        Dictionary with paths to generated files
    """
    visualizer = StrandVisualizer(autobreak_instance, override)
    return visualizer.generate_all_visualizations()


if __name__ == "__main__":
    # This would be called after autobreak processing is complete
    # create_strand_visualizations(my_autobreak_instance)
    pass