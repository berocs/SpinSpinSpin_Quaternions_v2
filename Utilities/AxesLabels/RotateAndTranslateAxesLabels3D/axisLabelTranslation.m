function                                                                     ...
axisLabelTranslation                                                         ...
    (                                                                        ...
      axesHandle,                                                            ...
      translationPoint,                                                      ...
      translationDirectionVector,                                            ...
      translationUnits                                                       ...
    )
%===============================================================================
%|
%|  FUNCTION:
%|    axisLabelTranslation
%|
%|  PURPOSE:
%|    This function moves axis labels to a proper distance from the axes.
%|
%|    It is called by function 'alignAxisLabels'.
%|
%|    The function 'alignAxisLabels' first rotates the x, y and z axis
%|    labels.
%|    Then this function is invoked to move those labels.
%|
%|    If you just want to move the labels but not rotate them,
%|    you have to specify the new position for the labels and
%|    the direction along which you move labels away from the
%|    axes.
%|
%|  INPUT(s):
%|
%|    axesHandle
%|      Handle of the axes.
%|
%|    translationPoint
%|      A 3-by-3 matrix,
%|      The first  column defines the position
%|      (in the original 3D space rather than in the 2D canvas space) 
%|      where you put the X axis label.
%|      The second column defines the position
%|      (in the original 3D space rather than in the 2D canvas space) 
%|      where you put the Y axis label.
%|      The third  column defines the position
%|      (in the original 3D space rather than in the 2D canvas space) 
%|      where you put the Z axis label.
%|
%|    translationDirectionVector
%|      The translation driection vector,
%|      a 2-by-3 matrix.
%|      The first  column defines the direction along which you move
%|      the x axis label away from the axis.
%|      The second column defines the direction along which you move
%|      the y axis label away from the axis.
%|      The third  column defines the direction along which you move
%|      the z axis label away from the axis.
%|      The vector should be normalized.
%|
%|    translationUnits
%|       'pixels',       Move labels away from axes in pixels.
%|       'characters',   Move labels away from axes in characters.
%|
%|  OUTPUT(s):
%|    None
%|
%===============================================================================


%===============================================================================
%        1         2         3         4         5         6         7         8
%2345678901234567890123456789012345678901234567890123456789012345678901234567890
%===============================================================================


%{------------------------------------------------------------------------------
   global AXISALIGN_TRANS_A;
   global AXISALIGN_TRANS_B;
%-------------------------------------------------------------------------------
   defaultValueA = 1.25;
   defaultValueB = 1.00;
%-------------------------------------------------------------------------------
   if( isempty( AXISALIGN_TRANS_A ) == true )
    %{--------------------------------------------------------------------------
    %  Global axis translation A is empty.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %  Assign default value.
    %---------------------------------------------------------------------------
       AXISALIGN_TRANS_A = defaultValueA;
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   if( isempty( AXISALIGN_TRANS_B ) == true )
    %{--------------------------------------------------------------------------
    %  Global axis translation B is empty.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %  Assign default value.
    %---------------------------------------------------------------------------
       AXISALIGN_TRANS_B = defaultValueB;
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   actualNumberInputArgs          = nargin;
   maximumExpectedNumberInputArgs = 4;
%-------------------------------------------------------------------------------
   if( actualNumberInputArgs < maximumExpectedNumberInputArgs )
    %{--------------------------------------------------------------------------
    %  There are fewer than the maximum number of input arguments.
    %  The translation units have not been specified.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %  Set the translation units to their default value (pixels).
    %---------------------------------------------------------------------------
       translationUnits = 'pixels';
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   matlabVersionInfo   = version( );
   matlabVersionNumber = str2double( matlabVersionInfo( 1 : 3 ) );
%-------------------------------------------------------------------------------
   xAxisLabelHandle = get( axesHandle, 'xlabel' );
   yAxisLabelHandle = get( axesHandle, 'ylabel' );
   zAxisLabelHandle = get( axesHandle, 'zlabel' );
%-------------------------------------------------------------------------------
   PIXEL_UNITS      = 1;
   CHARACTER_UNITS  = 2;
%-------------------------------------------------------------------------------
   if( strcmpi( translationUnits, 'pixels' ) == true )
    %{--------------------------------------------------------------------------
    %  Units are 'PIXELS'.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       translationUnits = PIXEL_UNITS;
    %---------------------------------------------------------------------------
       xLabelFontSize = get( xAxisLabelHandle, 'FontSize' )
    %---------------------------------------------------------------------------
    %  charPixelLength = 0.67 * xLabelFontSize;
    %  basePixelLength = 2.00 * xLabelFontSize;
    %---------------------------------------------------------------------------
       charPixelLength = 0.55; % 0.67 * xLabelFontSize;
       basePixelLength = 2.00; % 2.0   * xLabelFontSize;
    %}--------------------------------------------------------------------------
   elseif                                                                    ...
     ( strcmpi( translationUnits, 'character' ) == true )
    %{--------------------------------------------------------------------------
    %  Units are 'CHARACTERS'.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       translationUnits = CHARACTER_UNITS;
    %---------------------------------------------------------------------------
    %  Modify the method by which labels are moved.
    %---------------------------------------------------------------------------
    %  pos = pos +                                                           ...
    %        ( ( translationDirectionVector' ) .* aniso_ratio ) *            ...
    %        ( ( A * char_len ) + B );
    %---------------------------------------------------------------------------
       A           = AXISALIGN_TRANS_A;
       B           = AXISALIGN_TRANS_B;
       aniso_ratio = [ 2, 1 ];
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
%
%  Move X Label
%
%-------------------------------------------------------------------------------
   set(                                                                      ...
        xAxisLabelHandle,                                                    ...
        'HorizontalAlignment',          'Center',                            ...
        'VerticalAlignment',            'Middle',                            ...
        'Units',                        'data'                               ...
      );
   set(                                                                      ...
        xAxisLabelHandle,                                                    ...
        'Position',                     translationPoint( :, 1 )             ...
      );
   char_len_x = determineMaximumTickLabelLength                              ...
                         (                                                   ...
                           axesHandle,                                       ...
                           matlabVersionNumber,                              ...
                           'x'                                               ...
                         );
%-------------------------------------------------------------------------------
   if( translationUnits == CHARACTER_UNITS )
    %{--------------------------------------------------------------------------
    %  Units are CHARACTERS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            xAxisLabelHandle,                                                ...
            'Units',               'Characters'                              ...
          );
       xLabelPosition  = get(                                                ...
                              xAxisLabelHandle,  'Position'                  ...
                            );
       xLabelPosition( 1 : 2 ) =                                             ...
       xLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            (                                                ...
                              ( translationDirectionVector( :, 1 )' ) .*     ...
                              aniso_ratio                                    ...
                            ) *                                              ...
                            ( ( A * char_len_x ) + B )                       ...
                          );
    %}--------------------------------------------------------------------------
   else
    %{--------------------------------------------------------------------------
    %  Units are PIXELS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            xAxisLabelHandle,                                                ...
            'Units',               'Pixels'                                  ...
            );
       xLabelPosition  = get(                                                ...
                              xAxisLabelHandle,  'Position'                  ...
                            );
       xLabelPosition( 1 : 2 ) =                                             ...
       xLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            ( translationDirectionVector( :, 1 )' ) *        ...
                            (                                                ...
                              basePixelLength +                              ...
                              ( char_len_x * charPixelLength )               ...
                            )                                                ...
                          );
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   set(                                                                      ...
        xAxisLabelHandle,                                                    ...
        'Position',                     xLabelPosition                       ...
      );
%-------------------------------------------------------------------------------
%
%  Move Y Label
%
%-------------------------------------------------------------------------------
   set(                                                                      ...
        yAxisLabelHandle,                                                    ...
        'HorizontalAlignment',          'Center',                            ...
        'VerticalAlignment',            'Middle',                            ...
        'Units',                        'data'                               ...
      );
   set(                                                                      ...
        yAxisLabelHandle,                                                    ...
        'Position',                     translationPoint( :, 2 )             ...
      );
   char_len_y = determineMaximumTickLabelLength                              ...
                         (                                                   ...
                           axesHandle,                                       ...
                           matlabVersionNumber,                              ...
                           'y'                                               ...
                         );
%-------------------------------------------------------------------------------
   if( translationUnits == CHARACTER_UNITS )
    %{--------------------------------------------------------------------------
    %  Units are CHARACTERS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            yAxisLabelHandle,                                                ...
            'Units',               'Characters'                              ...
          );
       yLabelPosition  = get(                                                ...
                              yAxisLabelHandle,  'Position'                  ...
                            );
       yLabelPosition( 1 : 2 ) =                                             ...
       yLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            (                                                ...
                              ( translationDirectionVector( :, 2 )' ) .*     ...
                              aniso_ratio                                    ...
                            ) *                                              ...
                            ( ( A * char_len_y ) + B )                       ...
                          );
    %}--------------------------------------------------------------------------
   else
    %{--------------------------------------------------------------------------
    %  Units are PIXELS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            yAxisLabelHandle,                                                ...
            'Units',               'Pixels'                                  ...
            );
       yLabelPosition  = get(                                                ...
                              yAxisLabelHandle,  'Position'                  ...
                            );
       yLabelPosition( 1 : 2 ) =                                             ...
       yLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            ( translationDirectionVector( :, 2 )' ) *        ...
                            (                                                ...
                              basePixelLength +                              ...
                              ( char_len_y * charPixelLength )               ...
                            )                                                ...
                          );
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   set(                                                                      ...
        yAxisLabelHandle,                                                    ...
        'Position',                     yLabelPosition                       ...
      );
%-------------------------------------------------------------------------------
%
%  Move Z Label
%
%-------------------------------------------------------------------------------
   set(                                                                      ...
        zAxisLabelHandle,                                                    ...
        'HorizontalAlignment',          'Center',                            ...
        'VerticalAlignment',            'Bottom',                            ...
        'Units',                        'data'                               ...
      );
   set(                                                                      ...
        zAxisLabelHandle,                                                    ...
        'Position',                     translationPoint( :, 3 )             ...
      );
   char_len_z = determineMaximumTickLabelLength                              ...
                         (                                                   ...
                           axesHandle,                                       ...
                           matlabVersionNumber,                              ...
                           'z'                                               ...
                         );
%-------------------------------------------------------------------------------
   if( translationUnits == 2 )
    %{--------------------------------------------------------------------------
    %  Units are CHARACTERS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            zAxisLabelHandle,                                                ...
            'Units',               'Characters'                              ...
          );
       zLabelPosition  = get(                                                ...
                              zAxisLabelHandle,  'Position'                  ...
                            );
       zLabelPosition( 1 : 2 ) =                                             ...
       zLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            (                                                ...
                              ( translationDirectionVector( :, 3 )' ) .*     ...
                              aniso_ratio                                    ...
                            ) *                                              ...
                            ( ( A * char_len_z ) + B )                       ...
                          );
    %}--------------------------------------------------------------------------
   else
    %{--------------------------------------------------------------------------
    %  Units are PIXELS.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %---------------------------------------------------------------------------
       set(                                                                  ...
            zAxisLabelHandle,                                                ...
            'Units',               'Pixels'                                  ...
            );
       zLabelPosition  = get(                                                ...
                              zAxisLabelHandle,  'Position'                  ...
                            );
       zLabelPosition( 1 : 2 ) =                                             ...
       zLabelPosition( 1 : 2 ) +                                             ...
                          (                                                  ...
                            ( translationDirectionVector( :, 3 )' ) *        ...
                            (                                                ...
                              basePixelLength +                              ...
                              ( char_len_z * charPixelLength )               ...
                            )                                                ...
                          );
    %}--------------------------------------------------------------------------
   end;
%-------------------------------------------------------------------------------
   set(                                                                      ...
        zAxisLabelHandle,                                                    ...
        'Position',                     zLabelPosition                       ...
      );


%===============================================================================


    function                                                                 ...
    [                                                                        ...
      maximumTickLabelLength                                                 ...
    ] =                                                                      ...
    determineMaximumTickLabelLength                                          ...
             (                                                               ...
               axesHandle,                                                   ...
               matlabVersionNumber,                                          ...
               currentAxisName                                               ...
             )
    %{--------------------------------------------------------------------------
       [                                                                     ...
         tickLabel                                                           ...
       ] = get                                                               ...
            (                                                                ...
               axesHandle,                                                   ...
               [ currentAxisName,  'TickLabel' ]                             ...
            );
    %---------------------------------------------------------------------------
       if( matlabVersionNumber < 8.4 )
        %{----------------------------------------------------------------------
        %  Matlab version is previous to 8.4.
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %
        %-----------------------------------------------------------------------
           maximumTickLabelLength = size( tickLabel, 2 );
        %}----------------------------------------------------------------------
       else
        %{----------------------------------------------------------------------
        %  Matlab version is 8.4 or beyond.
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %  Must determine maximum tick length over all dimensions.
        %-----------------------------------------------------------------------
           maximumTickLabelLength = 0;
           numberDimensions       = length( tickLabel );
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           for( dimensionIndex = 1 : numberDimensions )
             %{-----------------------------------------------------------------
                currentTickLength = length( tickLabel{ dimensionIndex } );
             %------------------------------------------------------------------
                if( currentTickLength > maximumTickLabelLength )
                 %{-------------------------------------------------------------
                    maximumTickLabelLength = currentTickLength;
                 %}-------------------------------------------------------------
                end;
             %}-----------------------------------------------------------------
           end;
        %}----------------------------------------------------------------------
       end;
    %---------------------------------------------------------------------------
       return;
    %}--------------------------------------------------------------------------

%===============================================================================
   return;
%}------------------------------------------------------------------------------
