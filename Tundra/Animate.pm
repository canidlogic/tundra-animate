package Tundra::Animate;
use strict;

# Core dependencies
use Scalar::Util qw(looks_like_number);

# Non-core dependencies
use JSON::Tiny qw(decode_json);

=head1 NAME

Tundra::Animate - Interpret a JSON animation curves file.

=head1 SYNOPSIS

  use Tundra::Animate;
  
  # Parse the curves JSON files
  my $curves = Tundra::Animate->parse("curves.json");
  
  # Now generate all the frames
  for(my $f = 0; $f < $curves->frame_count; $f++) {
    
    # Generate all values for a specific frame
    my $ap = $curves->frame($f);
    
    # Now we can access curves for this frame
    my $param_value = $ap->{'my_param_name'};
    
    ...
  }

=head1 DESCRIPTION

The JSON curves file defines a set of custom-defined parameters from a
whole sequence of frames.  The file format supports linear interpolation
of parameter values across sequences of frames.

This module is able to parse the JSON curves file, and then generate
has references that are snapshots of each parameter value at requested
frame indices.

=head2 JSON format

The curves data file follows the JSON format.  The top-level entity in
the JSON must be a JSON object.  This object must have a property named
C<frame_count> that resolves to an integer storing the total number of
frames that will be animated, which must be at least one.

All other properties besides C<frame_count> are custom properties that
can be defined with any name and any meaning.  These other properties
are all the properties that will be generated within the hash reference
object for each specific frame.

The values of these custom properties must have a specific format.  The
value must be a JSON object that has properties C<vtype> and C<vcurve>.
The C<vtype> property is a string that names the data type of this
custom parameter.  It must be either C<"boolean">, C<"float">, or
C<"tuple">.

The C<"boolean"> type only supports C<true> and C<false> values.  The
C<"float"> type is a single, numeric floating-point value.  The
C<"tuple"> type is an array of two or more numeric floating-point
values, such as you might use to represent coordinates.

If C<vtype> is C<"tuple"> then there must also be a property C<vdim>
that is an integer defining how many numeric elements are in each array
value.  This must be at least two.  For example, you would use C<vdim>
with 2 for two-dimensional coordinates and C<vdim> with 3 for
three-dimensional coordinates.  You may also use V<vdim> values that are
greater than 3.

The C<vcurve> property defines the value of the parameter across all
frames.  The value of this property must be an array of one or more
I<shape arrays>.  Each shape array is an array of one or more I<sample
pairs>.  Each sample pair is an array of two values, the first of which
is a frame index (where zero is the first frame), and the second of
which is the value of this parameter at that frame.  The type of the
parameter value within the sample pair must match the C<vtype> declared
for this custom parameter.  For a C<vtype> of C<"boolean"> it must be
either true or false.  For a C<vtype> of C<"float"> it must be a numeric
value.  For a C<vtype> of C<"tuple"> it must be an array of numeric
values, the length of which matches the C<vdim> declared for this
parameter.

If a shape array contains only a single sample pair, then the shape
array defines the parameter value at the frame index given within that
sample pair.  If the shape array contains two or more sample pairs, then
the shape array defines the linear interpolation of the parameter value
starting at the frame index given within the first sample pair in the
shape, and proceeding up to and including the frame index given within
the last sample pair in the shape.  If there are more than two sample
pairs within the shape array, the values between the first and last
sample pairs represent intermediate parameter values.  In that case, the
curve looks like a sequence of connected line segments, where the first
sample pair is the start of the shape, the last sample pair is the end
of the shape, and all sample pairs in between are points that are
connected to each other with lines, resulting in a zig-zag shape after
linear interpolation is applied.

It is not possible to linearly interpolate a C<"boolean"> value.
Therefore, within in each shape array for a boolean parameter, each
sample pair must have the same boolean value or an error occurs during
parsing.

Within shape arrays, the frame indices of sample pairs must be in
strictly ascending order or an error will occur during parsing.  The
sequence of shape arrays must be continuous, such that if a shape array
is not the first shape array the frame index of the last sample pair in
the previous shape array must be exactly one less than the frame index
in the first sample pair of the current shape array, and if a shape
array is not the last shape array the frame index of the first sample
pair in the next shape array must be exactly one greater than the frame
index in the last sample pair of the current shape array.  Finally, the
frame index of the first sample pair of the first shape array must be
zero, and the frame index of the last sample pair of the last shape
array must be one less than the frame count.  Each parameter curve
therefore defines values for the full range of frames.

=head1 CONSTRUCTOR

=over 4

=item B<Tundra::Animate->parse(json_path)>

Construct a new Tundra Animate object instance by parsing the JSON file
at the given path.

=cut

sub parse {
  
  # Check parameter count
  ($#_ == 1) or die "Wrong number of parameters, stopped";
  
  # Get invocant and parameters
  my $invocant = shift;
  my $class = ref($invocant) || $invocant;
  
  my $json_path = shift;
  $json_path = "$json_path";
  
  # Check that JSON file exists
  (-f $json_path) or
    die "Can't find JSON curves file '$json_path', stopped";
  
  # Open the file in raw mode
  open(my $fh, "< :raw", $json_path) or
    die "Failed to open JSON curves file '$json_path', stopped";
  
  # Slurp the whole JSON file into memory
  my $js;
  {
    local $/;
    $js = readline($fh);
  }
  
  # Close the JSON file
  close($fh);
  
  # Parse the raw JSON
  $js = decode_json($js);
  
  # Top-level entity must be hash ref
  (ref($js) eq 'HASH') or
    die "Curves file must store JSON object, stopped";
  
  # Must have an unsigned integer frame_count that is greater than zero,
  # and convert it to integer
  (exists $js->{'frame_count'}) or
    die "Curves object must have frame_count property, stopped";
  
  (not ref($js->{'frame_count'})) or
    die "Curves frame_count property must be scalar, stopped";
  
  my $fcstr = $js->{'frame_count'};
  $fcstr = "$fcstr";
  ($fcstr =~ /^[1-9][0-9]*$/) or
    die "Curves frame_count must be non-zero unsigned integer, stopped";
  $js->{'frame_count'} = int($fcstr);
  
  # Now go through all other properties and verify their values
  for my $k (keys %$js) {
    
    # Skip verification for the special "frame_count" property that we
    # already verified
    if ($k eq 'frame_count') {
      next;
    }
    
    # Other properties must have object values
    my $v = $js->{$k};
    (ref($v) eq 'HASH') or
      die "Curves parameter must have JSON object value, stopped";
    
    # That object value must have vtype and vcurve properties
    for my $req_prop ("vtype", "vcurve") {
      (exists $v->{$req_prop}) or
        die "Curves parameter missing property '$req_prop', stopped";
    }
    
    # The vtype must be a scalar that is either "boolean", "float", or
    # "tuple"
    (not ref($v->{'vtype'})) or
      die "Curves vtype property must have scalar value, stopped";
    
    (($v->{'vtype'} eq 'boolean') or
      ($v->{'vtype'} eq 'float') or
      ($v->{'vtype'} eq 'tuple'))
        or die "Unrecognized curves vtype '$v->{'vtype'}', stopped";
    
    # If vtype is "tuple", then there must be a vdim property that has
    # an unsigned integer value that is two or greater; check that and
    # convert it to integer
    if ($v->{'vtype'} eq 'tuple') {
      (exists $v->{'vdim'}) or
        die "Curves tuple parameter requires vdim property, stopped";
      
      (not ref($v->{'vdim'})) or
        die "Curves vdim property must have scalar value, stopped";
      
      my $vdstr = $v->{'vdim'};
      $vdstr = "$vdstr";
      ($vdstr =~ /^[0-9]+$/) or
        die "Curves vdim property must be unsigned integer, stopped";
      
      $vdstr = int($vdstr);
      ($vdstr >= 2) or
        die "Curves vdim property must be at least two, stopped";
      
      $v->{'vdim'} = $vdstr;
    }
    
    # Check that the vcurve value is an array ref with at least one
    # element, that each element within this array is also an array ref
    # to an array of at least one element, and that each element within
    # that array is an array ref to an array of exactly two elements
    my $vc = $v->{'vcurve'};
    (ref($vc) eq 'ARRAY') or
      die "Curves vcurve property must be an array, stopped";
    
    (scalar(@$vc) >= 1) or
      die "Curves vcurve array may not be empty, stopped";
    
    for my $vca (@$vc) {
      (ref($vca) eq 'ARRAY') or
        die "Curves vcurve elements must be arrays, stopped";
      
      (scalar(@$vca) >= 1) or
        die "Curves vcurve elements must be non-empty, stopped";
      
      for my $vcb (@$vca) {
        (ref($vcb) eq 'ARRAY') or
          die "Curves vcurve elements must contain arrays, stopped";
        
        (scalar(@$vcb) == 2) or
          die "Curves vcurve element subarrays must be pairs, stopped";
      }
    }
    
    # Check that all sample values match the vtype and for tuple sample
    # values they have the correct dimension; for boolean sample values,
    # also make sure that each shape only has a single boolean value
    # across all its samples; also converts numerics to actual numbers
    # and converts booleans to 1 or 0
    for my $shape (@$vc) {
      my $first_sample = $shape->[0]->[1];
      
      for my $samp (@$shape) {
        my $samp_val = $samp->[1];
        
        if ($v->{'vtype'} eq 'tuple') {
          (ref($samp_val) eq 'ARRAY') or
            die "Sample type mismatch for $k, stopped";
          (scalar(@$samp_val) == $v->{'vdim'}) or
            die "Sample type misdimension for $k, stopped";
          for(my $ak = 0; $ak < scalar(@$samp_val); $ak++) {
            (not ref($samp_val->[$ak])) or
              die "Subsample type mismatch for $k, stopped";
            (looks_like_number($samp_val->[$ak])) or
              die "Subsample non-numeric for $k, stopped";
            $samp_val->[$ak] = $samp_val->[$ak] + 0.0;
          }
          
        } elsif ($v->{'vtype'} eq 'float') {
          (not ref($samp_val)) or
            die "Sample type mismatch for $k, stopped";
          (looks_like_number($samp_val)) or
            die "Sample non-numeric for $k, stopped";
          $samp->[1] = $samp_val + 0.0;
          
        } elsif ($v->{'vtype'} eq 'boolean') {
          (($samp_val == $JSON::Tiny::TRUE) or
              ($samp_val == $JSON::Tiny::FALSE)) or
            die "Sample type mismatch for $k, stopped";
          
          if ($samp_val == $JSON::Tiny::TRUE) {
            $samp->[1] = 1;
          } elsif ($samp_val == $JSON::Tiny::FALSE) {
            $samp->[1] = 0;
          }
          
        } else {
          # Shouldn't happen
          die "Unexpected, stopped";
        }
        
        if ($v->{'vtype'} eq 'boolean') {
          ($samp_val == $first_sample) or
            die "Boolean parameters can't be interpolated, stopped";
        }
      }
    }
    
    # Finally, check the frame indices within all the samples, and
    # convert them to numbers
    my $first_i = 0;
    for my $shape (@$vc) {
      # Check types, ranges, and convert all frame indices
      for my $samp (@$shape) {
        my $svi = $samp->[0];
        (not ref($svi)) or
          die "Curves sample index must be scalar, stopped";
        $svi = "$svi";
        ($svi =~ /^[0-9]+$/) or
          die "Curves sample index must be unsigned integer, stopped";
        $svi = int($svi);
        (($svi >= 0) and ($svi < $js->{'frame_count'})) or
          die "Curves sample index out of range, stopped";
        $samp->[0] = $svi;
      }
      
      # Check starting frame index
      ($shape->[0]->[0] == $first_i) or
        die "Curve for $k is missing frame $first_i, stopped";
      
      # Check strictly ascending
      for(my $si = 1; $si < scalar(@$shape); $si++) {
        ($shape->[$si]->[0] > $shape->[$si - 1]->[0]) or
          die "Curve for $k has frame indices out of order, stopped";
      }
      
      # Update first_i
      $first_i = $shape->[scalar(@$shape) - 1]->[0] + 1;
    }
    
    # Make sure we ended up exactly at the end of the frame sequence
    ($first_i == $js->{'frame_count'}) or
      die "Curve for $k doesn't match frame count, stopped";
  }
  
  # Define the object, store the JSON within it, bless it, and return
  my $self = { };
  bless($self, $class);
  
  $self->{'js'} = $js;
  return $self;
}

=back

=head1 INSTANCE METHODS

=over 4

=item B<object->frame_count()>

Return the total number of frames described by this curves object.

=cut

sub frame_count {
  ($#_ == 0) or die "Wrong number of parameters, stopped";
  my $self = shift;
  return $self->{'js'}->{'frame_count'};
}

=over 4

=item B<object->frame(i)>

Compute the value of all parameters at frame C<i> and return a hash
reference that contains all the computed parameter values.

The returned hash does not contain the special C<frame_count> parameter
but contains all other custom parameters that were defined in the JSON.
Boolean values will map to 1 for true or 0 for false, float values will
map to the actual numeric value, and tuple values will map to an array
reference storing all of the numeric values.

=cut

sub frame {
  # Check parameter count
  ($#_ == 1) or die "Wrong number of parameters, stopped";
  
  # Get parameters
  my $self = shift;
  my $frame_i = shift;
  
  $frame_i = int($frame_i);
  
  # Check frame index
  (($frame_i >= 0) and ($frame_i < $self->{'js'}->{'frame_count'})) or
    die "Frame index out of range, stopped";
  
  # Declare hash that will hold computed values
  my %vh;
  
  # Go through all keys
  for my $k (keys %{$self->{'js'}}) {
    
    # Skip the special "frame_count" key
    if ($k eq 'frame_count') {
      next;
    }
    
    # Within the curves property, find the shape that contains this
    # frame index
    my $found_frame = 0;
    my $frame_shape;
    for my $shape (@{$self->{'js'}->{$k}->{'vcurve'}}) {
      if (($shape->[0]->[0] <= $frame_i) and
            ($shape->[scalar(@$shape) - 1]->[0] >= $frame_i)) {
        $found_frame = 1;
        $frame_shape = $shape;
        last;
      }
    }
    ($found_frame) or die "Unexpected, stopped at";
    
    # If this property is boolean, we don't interpolate; instead, just
    # take the value from the start of the shape and then continue to
    # the next key
    if ($self->{'js'}->{$k}->{'vtype'} eq 'boolean') {
      $vh{$k} = $frame_shape->[0]->[1];
      next;
    }
    
    # Not boolean, so now we need to find the last sample in the shape
    # that has a frame index less than or equal to the requested frame
    my $j;
    for($j = scalar(@$frame_shape) - 1; $j >= 0; $j--) {
      if ($frame_shape->[$j]->[0] <= $frame_i) {
        last;
      }
    }
    ($j >= 0) or die "Unexpected, stopped at";
    
    # If the requested frame is equal to the frame in the selected
    # sample, then add the value in the selected sample and continue
    # to next property without interpolation required
    if ($frame_shape->[$j]->[0] == $frame_i) {
      $vh{$k} = $frame_shape->[$j]->[1];
      next;
    }
    
    # If we got here, linear interpolation is required; j should never
    # point to the last element of a shape in this case
    ($j < scalar(@$frame_shape) - 1) or die "Unexpected, stopped at";
    
    # Compute the progress between this sample and the next so we know
    # how much to interpolate
    my $p = ($frame_i - $frame_shape->[$j]->[0]) /
              ($frame_shape->[$j + 1]->[0] - $frame_shape->[$j]->[0]);
    
    # Compute the interpolated value and add it to the results; we
    # already handled booleans earlier so we don't need to worry about
    # that here
    my $a = $frame_shape->[$j]->[1];
    my $b = $frame_shape->[$j + 1]->[1];
    my $ipv;
    
    if ($self->{'js'}->{$k}->{'vtype'} eq 'float') {
      $ipv = ($a * (1.0 - $p)) + ($b * $p);
      
    } elsif ($self->{'js'}->{$k}->{'vtype'} eq 'tuple') {
      my @aa = @$a;
      my @bb = @$b;
      $ipv = [];
      
      for(my $z = 0; $z <= $#aa; $z++) {
        my $cv = ($aa[$z] * (1.0 - $p)) + ($bb[$z] * $p);
        push @$ipv, ($cv);
      }
      
    } else {
      die "Unexpected, stopped at";
    }
    
    $vh{$k} = $ipv;
  }
  
  # Return a reference to the hash we computed
  return \%vh;
}

=back

=head1 AUTHOR

Noah Johnson, C<noah.johnson@loupmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2022 Multimedia Data Technology Inc.

MIT License:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=cut

# End with something that evaluates to true
#
1;
