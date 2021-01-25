<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Structured grid generation from NekMesh wrapper</description>
    <executable python="true">StructuredGrid.py</executable>
    <parameters> 5 6 0.0 1.0  2.0 3.0  2 Quad  output.xml </parameters>
    <files />
    <metrics>
        <metric type="regex" id="1">
            <regex>^.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
