<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Do 3D smoothing </description>
    <executable>FieldConvert</executable>
    <parameters> -m C0Projection:helmsmoothing=0.6 Hex_channel_C0helmsmoothing.xml Hex_channel_C0helmsmoothing.fld C0Hel.plt -e -f -v </parameters>
    <files>
        <file description="Session File">Hex_channel_C0helmsmoothing.xml</file>
        <file description="Session File">Hex_channel_C0helmsmoothing.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.57735</value>
            <value variable="y" tolerance="1e-4">0.57735</value>
            <value variable="z" tolerance="1e-4">0.57735</value>
            <value variable="u" tolerance="1e-6">0.635792</value>
            <value variable="v" tolerance="1e-6">0.635792</value>
            <value variable="w" tolerance="1e-6">0.635792</value>
            <value variable="p" tolerance="1e-6">1.10303</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="x" tolerance="1e-4">1.0</value>
            <value variable="y" tolerance="1e-4">1.0</value>
            <value variable="z" tolerance="1e-4">1.0</value>
            <value variable="u" tolerance="1e-6">0.739606</value>
            <value variable="v" tolerance="1e-6">0.739606</value>
            <value variable="w" tolerance="1e-6">0.739606</value>
            <value variable="p" tolerance="1e-6">2.21882</value>
        </metric>
    </metrics>
</test>

