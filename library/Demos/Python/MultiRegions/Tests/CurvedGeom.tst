<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Integral in standard 2D quad region </description>
    <executable python="true"> CurvedGeom.py </executable>
    <parameters></parameters>
    <metrics>
        <metric type="regex" id="1">
            <regex>L infinity error: ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)$</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-12">0.000000e+00</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
