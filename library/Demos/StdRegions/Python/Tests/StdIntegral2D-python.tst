<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Std Projection in 1D  </description>
    <executable python="true"> StdProject1D.py </executable>
    <parameters></parameters>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*error = ([-+]?\d*\.?\d+[eE]-+]?\d+)?$</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-12">4.44089e-16</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
