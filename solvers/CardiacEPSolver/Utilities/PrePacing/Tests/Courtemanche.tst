<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Courtemanche Cell model</description>
    <executable>PrePacing</executable>
    <parameters>Courtemanche.xml</parameters>
    <files>
        <file description="Session File">Courtemanche.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>
                ^#\s([\w]*)\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
            </regex>
            <matches>
                <match>
                    <field>u</field>
                    <field tolerance="1e-06">-7.51776</field>
                </match>
                <match>
                    <field>m</field>
                    <field tolerance="1e-06">0.987892</field>
                </match>
                <match>
                    <field>h</field>
                    <field tolerance="1e-06">2.03226e-178</field>
                </match>
                <match>
                    <field>j</field>
                    <field tolerance="1e-06">3.67039e-12</field>
                </match>
                <match>
                    <field>o_a</field>
                    <field tolerance="1e-06">0.677556</field>
                </match>
                <match>
                    <field>o_i</field>
                    <field tolerance="1e-06">0.00173797</field>
                </match>
                <match>
                    <field>u_a</field>
                    <field tolerance="1e-06">0.915341</field>
                </match>
                <match>
                    <field>u_i</field>
                    <field tolerance="1e-06">0.993285</field>
                </match>
                <match>
                    <field>x_r</field>
                    <field tolerance="1e-06">0.213132</field>
                </match>
                <match>
                    <field>x_s</field>
                    <field tolerance="1e-06">0.0846056</field>
                </match>
                <match>
                    <field>d</field>
                    <field tolerance="1e-06">0.580487</field>
                </match>
                <match>
                    <field>f</field>
                    <field tolerance="1e-06">0.67841</field>
                </match>
                <match>
                    <field>f_Ca</field>
                    <field tolerance="1e-06">0.350083</field>
                </match>
                <match>
                    <field>U</field>
                    <field tolerance="1e-06">0.999994</field>
                </match>
                <match>
                    <field>V</field>
                    <field tolerance="1e-06">2.5819e-11</field>
                </match>
                <match>
                    <field>W</field>
                    <field tolerance="1e-06">0.94222</field>
                </match>
                <match>
                    <field>Na_i</field>
                    <field tolerance="1e-06">11.1712</field>
                </match>
                <match>
                    <field>Ca_i</field>
                    <field tolerance="1e-06">0.000645922</field>
                </match>
                <match>
                    <field>K_i</field>
                    <field tolerance="1e-06">138.991</field>
                </match>
                <match>
                    <field>Ca_rel</field>
                    <field tolerance="1e-06">0.218687</field>
                </match>
                <match>
                    <field>Ca_up</field>
                    <field tolerance="1e-06">1.58454</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>




