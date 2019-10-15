<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) high-order tri square in circle with surface extraction test</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m extract:surf=2 extract_square_circ.msh "extract_square_circ-out.xml:xml:test:testcond=abs(x^2+y^2-0.0625)&lt;1e-4"</parameters>
    <files>
        <file description="Input File">extract_square_circ.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
